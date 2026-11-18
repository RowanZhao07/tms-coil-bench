function coils = build_coils(nPerTurn, rootDir, writeManifest)
%BUILD_COILS Fourteen documented, reproducible filament coil models in SI.
% segments = [ax ay az bx by bz currentWeight]; the operating current is
% supplied separately by the field solver. Negative weights reverse current.
% Dimensions given as inequalities in the source drawings are bounds, not exact
% centreline radii. Missing dimensions below are explicit modelling choices.
% All generated loops close. The imported measured H1 retains source gaps.
% All endpoints have z <= 0. Observation depth is positive z from this plane.
if nargin < 1 || isempty(nPerTurn), nPerTurn = 128; end
if nargin < 2 || isempty(rootDir), rootDir = fileparts(fileparts(mfilename('fullpath'))); end
if nargin < 3 || isempty(writeManifest), writeManifest = false; end
validateattributes(nPerTurn, {'numeric'}, {'scalar','integer','>=',32});
coils = struct('id',{},'name_en',{},'name_cn',{},'source',{},'note',{},'segments',{});

% 1. Circular: 9 turns in each of two layers. Tube OD=4 mm is used to
% interpret the clear bore/outer envelope; a filament itself has no thickness.
S = rings(linspace(17,55,9), [0 -4.5], 0, 0, nPerTurn);
coils(end+1) = item('circular','Circular','圆形线圈','Coil geometry drawings PDF p.1', ...
    'Two layers x 9 turns. Centreline radii 17:4.75:55 mm; layer pitch 4.5 mm. With assumed 4-mm conductor OD: clear bore 30 mm, OD 114<115 mm, thickness 8.5<10 mm. Feed and spiral interconnections omitted; closed concentric-turn equivalent.',S);

% 2. H1: actual public X-ray reconstructed geometry. Proper rotation Rx(pi).
h = load(fullfile(rootDir,'reference_models','h1_segments_si.mat'));
R = diag([1 -1 -1]); a=h.start_m*R.'; b=h.end_m*R.';
allp=[a;b]; shift=[-(min(allp(:,1))+max(allp(:,1)))/2, ...
    -(min(allp(:,2))+max(allp(:,2)))/2, -max(allp(:,3))];
a=a+shift; b=b+shift;
S=[a b ones(size(a,1),1)];
coils(end+1) = item('h1','H1 (X-ray reconstruction)','H1线圈', ...
    'SimNIBS Brainsway_H1.tcd v1.0.o; commit fde509f6f8de38a4b03947944cf7ee93b8bebabe; Worbs et al. 2025 doi:10.1016/j.brs.2025.05.136', ...
    sprintf(['14895 source line elements, 18 components, default deformations 0 deg. Proper rotation Rx(pi)=diag(1,-1,-1), translation [%.9g %.9g %.9g] m. ', ...
    'Endpoint bounding box centred in x/y; uppermost endpoint z=0. This is a common air-field placement, NOT an anatomical fit. Original total dl=[-3.345800,2.643905,0.463064] mm: small source gaps retained, no fictitious return wire added. Source sampling is fixed and is not refined by nPerTurn.'],shift),S);

% 3. Simplified H1: bilateral scalp-following active arcs and distant returns.
% No dimensions are supplied for this schematic in the source drawings.
S=zeros(0,7);
for side=[-1 1]
    for k=0:5
        th=linspace(0,side*75*pi/180,max(16,round(nPerTurn/3))).';
        radius=90+2*k;
        active=[radius*sin(th), -25+5*k+zeros(size(th)), radius*(cos(th)-1)];
        endpoint=active(end,:);
        P=[active; endpoint(1),55+5*k,-140-3*k; 0,55+5*k,-140-3*k; active(1,:)];
        S=[S; polysegments(P,side)]; %#ok<AGROW>
    end
end
coils(end+1) = item('h1_simplified','Simplified H1 (schematic)','简化H1线圈','Coil geometry drawings PDF p.2 schematic; dimensions unspecified', ...
    'Assumed bilateral 6-loop design. Active arc radii 90..100 mm, 75 deg per side, y=-25..0 mm. Return legs y=55..80 mm, z=-140..-155 mm. Closed 3D loops, opposed left/right circulation reinforces the active arcs. This is an engineering schematic, not a validated Brainsway H1 replica.',S);

% 4. Alternative H1: full transverse active arches with raised rectangular
% return conductors, following the qualitative silhouette on PDF p.3.
S=zeros(0,7);
for k=0:7
    th=linspace(-72*pi/180,72*pi/180,max(24,round(nPerTurn/2))).';
    radius=84+2*k;
    active=[radius*sin(th), -28+4*k+zeros(size(th)), radius*(cos(th)-1)];
    P=[active; active(end,1),60+4*k,-150-3*k; active(1,1),60+4*k,-150-3*k; active(1,:)];
    S=[S; polysegments(P,1)]; %#ok<AGROW>
end
coils(end+1) = item('h1_alternative','Alternative H1 (schematic)','另类H1线圈','Coil geometry drawings PDF p.3 schematic; dimensions unspecified', ...
    'Assumed 8 closed transverse loops. Active arch radii 84..98 mm, +/-72 deg, y=-28..0 mm; return y=60..88 mm and z=-150..-171 mm. Qualitative arch/rectangular-return interpretation only, not an identified commercial H1 geometry.',S);

% 5. The source drawings repeat "children figure eight" with inconsistent sizes.
% Use the explicitly annotated compact four-layer drawing (p.9) as primary.
S=pairrings([17 21.5 26 30.5], -(0:3)*4.5, 35, nPerTurn);
coils(end+1) = item('pediatric_figure8','Pediatric figure-of-eight','儿童8字线圈','Coil geometry drawings PDF p.9 (primary); conflicting larger drawing p.3-4', ...
    '4 turns per lobe per layer, 4 layers, centres +/-35 mm. Radii [17,21.5,26,30.5] mm and layer pitch 4.5 mm assume 4-mm OD: bore30, lobe OD65<70, width135<136, thickness17.5<21 mm. Opposite lobe current. Earlier p.3-4 drawing instead has centre distance92 and OD90/bore25; it is a separate larger variant, not silently merged.',S);

% 6. "Elliptical" drawing is a racetrack, not a mathematical ellipse.
S=zeros(0,7);
for z=[0 -4.5]
    for r=17:4:45
        th1=linspace(0,pi,max(16,round(nPerTurn/2))).';
        th2=linspace(pi,2*pi,max(16,round(nPerTurn/2))).';
        P=[r*cos(th1),40+r*sin(th1),z+zeros(size(th1)); ...
           r*cos(th2),-40+r*sin(th2),z+zeros(size(th2))];
        S=[S; polysegments(P,1)]; %#ok<AGROW>
    end
end
coils(end+1) = item('elliptical','Elliptical / racetrack','椭圆线圈','Coil geometry drawings PDF p.5', ...
    'Racetrack interpretation: two straight sections of 80 mm plus semicircular ends. Assumed 8 turns x 2 layers, end radii17:4:45 mm, pitch4.5 mm between layers. 4-mm OD gives inner30x110 and outer94x174 mm, within <100x175 mm; thickness8.5<12 mm. Turn count is a modelling choice because the drawing does not state it.',S);

% 7. Double cone: exact 95-degree included angle via symmetric 42.5-degree
% rigid wing rotations about the central y axis. Leads are omitted.
S=pairrings(linspace(31.1,50,5),[0 -4.5],52,nPerTurn);
S=wing_rotate(S,42.5);
coils(end+1) = item('double_cone','Double cone','双锥线圈','Coil geometry drawings PDF p.6', ...
    '5 turns x 2 layers per wing; radii31.1..50 mm, nominal unbent centres +/-52 mm. Rigid wings tilted42.5 deg from plane -> included95 deg. 4-mm OD gives clear bore58.2>58 mm and OD104<107 mm. Drawing also marks overall190 mm without a unique hinge/centre definition; centre distance104 mm is assumed, so this model does not claim to reproduce that overall dimension.',S);

% 8. Single-layer 18-turn pancake (outer bound, not required OD=170 mm).
S=rings(7:4:75,0,0,0,nPerTurn);
coils(end+1) = item('pancake','Bottom pancake','盘底线圈','Coil geometry drawings PDF p.7', ...
    '18 turns x 1 layer; centreline radii7:4:75 mm. 4-mm tube OD gives clear bore10 mm, outer radius77 mm (OD154<170 mm), matching the drawing handwritten outer-radius interpretation. Source nominal inductance20+/-3 uH is not imposed or verified by this B-field model.',S);

% 9. Fractional source turns represented explicitly as an ampere-turn
% equivalent. It is not a geometrically faithful half-turn and return lead.
S=zeros(0,7);
for z=-(0:3)*4.5
    for k=0:7
        w=1; if k==7,w=0.5;end
        P=circle(9.5+4*k,0,0,z,nPerTurn);
        S=[S; polysegments(P,w)]; %#ok<AGROW>
    end
end
coils(end+1) = item('pediatric_circular','Pediatric circular','儿童圆形线圈','Coil geometry drawings PDF p.8', ...
    '4 layers x 7.5 effective turns: 7 closed unit-current rings plus outermost closed ring weighted0.5 in each layer (30 ampere-turns at1A). Radii9.5:4:37.5 mm; layer pitch4.5 mm. 4-mm OD gives bore15, OD79<80, thickness17.5<20 mm. The 0.5 weighting is an explicit ampere-turn surrogate, not the real 7.5-turn spiral/return routing.',S);

% 10. T refers to the out-of-plane lead/handle in p.10, not a third lobe.
S=pairrings(linspace(14.5,44.5,9),[0 -4.5],46,nPerTurn);
coils(end+1) = item('t_coil','T-type (winding equivalent)','T型线圈','Coil geometry drawings PDF p.10; p.9 repeats a different drawing under the heading', ...
    'Two planar lobes, opposed current, centres +/-46 mm, assumed9 turns x2 layers, radii14.5..44.5 mm; radial pitch3.75 mm and layer pitch4.5 mm. Conductor OD is not labelled on p.10: assumed3 mm gives bore26, totalwidth184<185 mm and thickness7.5<12 mm, without radial overlap. Actual T-shaped out-of-plane leads shown on p.10 are omitted because complete supply/return routing is missing. This closed-turn winding equivalent cannot distinguish handle-only variants.',S);

% 11. Two-stage butterfly folding: central included angle160 deg and an
% additional43-deg fold (supplement of137 deg) beyond each lobe centre.
S=pairrings(linspace(31.1,50,5),[0 -4.5],53,nPerTurn);
S=butterfly_fold(S,53,10,43);
coils(end+1) = item('butterfly','Butterfly (two-stage folded)','蝶形线圈','Coil geometry drawings PDF p.11', ...
    '5 turns x2 layers per lobe; radii31.1..50 mm, nominal unbent centre distance106 mm. Central wings10 deg (included160 deg); outer halves add43 deg (137-deg internal bend), hinge at lobe centre53 mm. This is an explicitly chosen interpretation of the compound-bend sketch. 4-mm OD gives bore58.2>58 and lobeOD104<105 mm; finite casing/clearances are not modelled.',S);

% 12. Three symmetric figure-eight pairs fan upwards, as PDF p.12 depicts.
S=zeros(0,7); tt=linspace(0,2*pi,nPerTurn+1).';
for angle=[0 35 65]
    for side=[-1 1]
        a=44; b=28;
        for k=0:5
            aa=a-4*k; bb=b-3*k;
            P=[a+aa*cos(tt),bb*sin(tt),zeros(size(tt))];
            rot=angle*pi/180; M=[cos(rot) -sin(rot);sin(rot) cos(rot)];
            P(:,1:2)=P(:,1:2)*M.';
            P(:,1)=side*P(:,1);
            % Mirroring already reverses circulation; use equal weights.
            S=[S; polysegments(P,1)]; %#ok<AGROW>
        end
    end
end
coils(end+1) = item('figure8_array','Figure-of-eight array','8字阵列线圈','Coil geometry drawings PDF p.12 schematic; no dimensions/angles/current phases specified', ...
    'Assumed 3 mirror-symmetric pairs fanning at0,35,65 deg. Each lobe6 nested ellipses, fixed centre44 mm from origin, semi-axes44-4k and28-3k mm, k0..5. Mirrors reverse circulation. All three pair channels are synchronous equal amplitude; this is one explicit array drive, not a uniquely prescribed optimal array. Crossing filament intersections are idealized; no manufacturing clearance claimed.',S);

% 13. Only the outer wing portions are folded, not both whole lobes.
S=pairrings(linspace(20,44,6),0,46,nPerTurn);
S=outer_fold(S,56,45,45);
coils(end+1) = item('folded_figure8','Folded figure-of-eight','折叠8字线圈','Coil geometry drawings PDF p.13-15 qualitative fold equations', ...
    'Assumed6 turns/lobe, radii20..44 mm, centres +/-46 mm; outer hinges x=+/-56 mm (dL=dR=10 mm from respective centres), betaL=betaR=45 deg away from observation side. Inner parts remain planar. Opposite lobe currents as p.14. No metric dimensions or fold angles are specified in the document.',S);

% 14. Two outer parts of the same circular winding fold about parallel chords.
S=rings(linspace(17,55,9),[0 -4.5],0,0,nPerTurn);
S=outer_fold(S,20,30,30);
coils(end+1) = item('folded_circular','Folded circular','折叠圆形线圈','Coil geometry drawings PDF p.16 schematic; no metric dimensions', ...
    'Uses circular reference radii17..55 mm,9 turns x2 layers; folds beyond x=+/-20 mm at beta1=beta2=30 deg away from observation side. Inner portion remains planar. Radius, hinge distances and fold angles are explicit assumptions; source only supplies symbolic d1,d2,beta1,beta2.',S);

% Normalise each geometry to the same empty observation half-space. This
% changes only placement, never size. Retain current weights and handedness.
for k=1:numel(coils)
    S=coils(k).segments;
    zmax=max([S(:,3);S(:,6)]); S(:,[3 6])=S(:,[3 6])-zmax;
    coils(k).segments=S;
    if any(~isfinite(S),'all') || any(vecnorm(S(:,4:6)-S(:,1:3),2,2)<1e-14)
        error('build_coils:InvalidSegment','Invalid segment in %s.',coils(k).id);
    end
end

% Save inspectable metadata without repeating tens of thousands of vertices.
manifest=rmfield(coils,'segments');
for k=1:numel(coils)
    S=coils(k).segments; pts=[S(:,1:3);S(:,4:6)];
    manifest(k).segment_count=size(S,1);
    manifest(k).bbox_m=[min(pts,[],1);max(pts,[],1)];
    manifest(k).weighted_length_m=sum(vecnorm(S(:,4:6)-S(:,1:3),2,2).*abs(S(:,7)));
    manifest(k).weighted_closure_vector_m=sum((S(:,4:6)-S(:,1:3)).*S(:,7),1);
end
if writeManifest
    fid=fopen(fullfile(rootDir,'geometry_manifest.json'),'w','n','UTF-8');
    if fid>=0, fprintf(fid,'%s',jsonencode(manifest,'PrettyPrint',true)); fclose(fid); end
end
end

function c=item(id,en,cn,source,note,S)
c=struct('id',id,'name_en',en,'name_cn',cn,'source',source,'note',note,'segments',S);
end

function P=circle(r,cx,cy,z,n)
t=linspace(0,2*pi,n+1).'; P=[cx+r*cos(t),cy+r*sin(t),z+zeros(size(t))];
end

function S=polysegments(P,w)
% Input millimetres, output metres; always explicitly close the path.
if norm(P(1,:)-P(end,:))>1e-9, P=[P;P(1,:)]; else, P(end,:)=P(1,:); end
a=P(1:end-1,:)*1e-3; b=P(2:end,:)*1e-3;
keep=vecnorm(b-a,2,2)>1e-14;
S=[a(keep,:),b(keep,:),w*ones(nnz(keep),1)];
end

function S=rings(radii,zs,cx,cy,n)
S=zeros(0,7);
for z=zs
    for r=radii, S=[S;polysegments(circle(r,cx,cy,z,n),1)]; end %#ok<AGROW>
end
end

function S=pairrings(radii,zs,c,n)
left=rings(radii,zs,-c,0,n); left(:,7)=-1;
right=rings(radii,zs,c,0,n); S=[left;right];
end

function S=wing_rotate(S,deg)
% Rotate each whole lobe by the sign of its centre. Never reflect geometry.
for cols={1:3,4:6}
    q=cols{1}; P=S(:,q); side=sign((S(:,1)+S(:,4))/2);
    ang=side*deg*pi/180; x=P(:,1); z=P(:,3);
    P(:,1)=cos(ang).*x+sin(ang).*z;
    P(:,3)=-sin(ang).*x+cos(ang).*z; S(:,q)=P;
end
end

function S=outer_fold(S,hinge_mm,left_deg,right_deg)
h=hinge_mm*1e-3;
for cols={1:3,4:6}
    q=cols{1}; P=S(:,q); x=P(:,1); d=max(abs(x)-h,0);
    deg=right_deg+zeros(size(x)); deg(x<0)=left_deg;
    a=deg*pi/180;
    P(:,1)=sign(x).*(min(abs(x),h)+d.*cos(a));
    P(:,3)=P(:,3)-d.*sin(a); S(:,q)=P;
end
end

function S=butterfly_fold(S,hinge_mm,inner_deg,extra_deg)
h=hinge_mm*1e-3; a=inner_deg*pi/180; b=(inner_deg+extra_deg)*pi/180;
for cols={1:3,4:6}
    q=cols{1}; P=S(:,q); x=P(:,1); d=abs(x); u=min(d,h); v=max(d-h,0);
    P(:,1)=sign(x).*(u*cos(a)+v*cos(b));
    P(:,3)=P(:,3)-u*sin(a)-v*sin(b); S(:,q)=P;
end
end
