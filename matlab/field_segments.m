function B = field_segments(segments, observers, current_A, method, useGPU)
%FIELD_SEGMENTS Magnetic flux density [T] from directed straight wires [m].
% segments = [ax ay az bx by bz signed_weight]; observers = N-by-3.
% exact: closed-form integral on each straight segment, not a dipole model.
% midpoint: independent midpoint quadrature on the same polygonal geometry.
% Chunked double precision keeps workspace well below an 8 GB GPU budget.
if nargin<4, method='exact'; end
if nargin<5, useGPU=false; end
assert(size(segments,2)==7 && size(observers,2)==3);
assert(all(isfinite(segments(:))) && all(isfinite(observers(:))));
mu0=4*pi*1e-7; % conventional engineering value, recorded in run metadata
a=segments(:,1:3); dl=segments(:,4:6)-a;
len=sqrt(sum(dl.^2,2));
keep=len>1e-14 & segments(:,7)~=0;
a=a(keep,:); dl=dl(keep,:); len=len(keep); w=segments(keep,7);
u=dl./len;
B=zeros(size(observers));
pointChunk=256; segmentChunk=512;
for p0=1:pointChunk:size(observers,1)
    ix=p0:min(p0+pointChunk-1,size(observers,1));
    p=observers(ix,:); accum=zeros(numel(ix),3);
    if useGPU, p=gpuArray(p); accum=gpuArray(accum); end
    for s0=1:segmentChunk:size(a,1)
        j=s0:min(s0+segmentChunk-1,size(a,1));
        aa=a(j,:); uu=u(j,:); ll=len(j)'; ww=w(j)'; dd=dl(j,:);
        if useGPU
            aa=gpuArray(aa); uu=gpuArray(uu); ll=gpuArray(ll);
            ww=gpuArray(ww); dd=gpuArray(dd);
        end
        rx=p(:,1)-aa(:,1)'; ry=p(:,2)-aa(:,2)'; rz=p(:,3)-aa(:,3)';
        if strcmp(method,'midpoint')
            rx=rx-dd(:,1)'/2; ry=ry-dd(:,2)'/2; rz=rz-dd(:,3)'/2;
            r2=rx.^2+ry.^2+rz.^2;
            if any(r2(:)<1e-24), error('Observer coincides with wire midpoint.'); end
            fac=ww./(r2.*sqrt(r2));
            cx=dd(:,2)'.*rz-dd(:,3)'.*ry;
            cy=dd(:,3)'.*rx-dd(:,1)'.*rz;
            cz=dd(:,1)'.*ry-dd(:,2)'.*rx;
        elseif strcmp(method,'exact')
            proj=rx.*uu(:,1)'+ry.*uu(:,2)'+rz.*uu(:,3)';
            qx=rx-proj.*uu(:,1)'; qy=ry-proj.*uu(:,2)'; qz=rz-proj.*uu(:,3)';
            rho2=qx.^2+qy.^2+qz.^2;
            onLine=rho2<1e-24;
            onWire=onLine & proj>=0 & proj<=ll;
            if any(onWire(:)), error('Observer lies on filament; finite-wire model required.'); end
            rho2(onLine)=1; % collinear outside the segment contributes exactly zero
            fac=ww./rho2.*(proj./sqrt(rho2+proj.^2)-(proj-ll)./sqrt(rho2+(proj-ll).^2));
            fac(onLine)=0;
            cx=uu(:,2)'.*qz-uu(:,3)'.*qy;
            cy=uu(:,3)'.*qx-uu(:,1)'.*qz;
            cz=uu(:,1)'.*qy-uu(:,2)'.*qx;
        else
            error('Unknown integration method: %s',method);
        end
        accum=accum+[sum(cx.*fac,2),sum(cy.*fac,2),sum(cz.*fac,2)];
    end
    if useGPU, accum=gather(accum); end
    B(ix,:)=mu0*current_A/(4*pi)*accum;
end
end
