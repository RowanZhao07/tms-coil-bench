// Axisymmetric magnetostatic transformed vector potential u=A_phi/r, SI.
// Exact 3D axisymmetry reduction: -div(r^3 grad u) = mu0 r^2 J_phi.
// Axis is physical z (Gmsh Y), radius r is Gmsh X.
// B_r=-r du/dz; B_z=2u+r du/dr. This avoids dividing by r on axis.
Include "case_parameters.pro";
Group {
  Air = Region[1]; Coil = Region[2];
  Domain = Region[{1,2}]; Boundary = Region[3];
}
Function {
  weight[Domain] = X[]^3;
  source[Coil] = (4.e-7*Pi)*X[]^2*CurrentDensity;
}
Constraint {
  { Name Dirichlet; Case { { Region Boundary; Value 0.; } } }
}
Jacobian { { Name J; Case { { Region All; Jacobian Vol; } } } }
Integration {
  { Name I; Case { { Type Gauss;
    Case { { GeoElement Triangle; NumberOfPoints 7; } }
  } } }
}
FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef un; Function BF_Node;
        Support Domain; Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef un; EntityType NodesOf; NameOfConstraint Dirichlet; }
    }
  }
}
Formulation {
  { Name Magnetostatics; Type FemEquation;
    Quantity { { Name u; Type Local; NameOfSpace Hgrad; } }
    Equation {
      Integral { [ weight[] * Dof{d u}, {d u} ]; In Domain; Jacobian J; Integration I; }
      Integral { [ -source[], {u} ]; In Coil; Jacobian J; Integration I; }
    }
  }
}
Resolution {
  { Name SolveA;
    System { { Name A; NameOfFormulation Magnetostatics; } }
    Operation { Generate[A]; Solve[A]; SaveSolution[A]; }
  }
}
PostProcessing {
  { Name Fields; NameOfFormulation Magnetostatics;
    Quantity {
      { Name b; Value { Term { [ Vector[-X[]*CompY[{d u}],2*{u}+X[]*CompX[{d u}],0] ]; In Domain; Jacobian J; } } }
      { Name a; Value { Term { [ X[]*{u} ]; In Domain; Jacobian J; } } }
    }
  }
}
PostOperation {
  { Name Export; NameOfPostProcessing Fields;
    Operation {
      Print[b, OnLine {{0,0.01,0}{0,0.05,0}} {4}, File "axis_b.txt", Format Table];
      Print[b, OnElementsOf Domain, File "field_b.pos"];
      Print[a, OnElementsOf Domain, File "field_a.pos"];
    }
  }
}
