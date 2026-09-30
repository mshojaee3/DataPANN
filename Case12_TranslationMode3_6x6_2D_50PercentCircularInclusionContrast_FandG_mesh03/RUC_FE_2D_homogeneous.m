function [T, diagOut] = RUC_FE_2D_homogeneous(H, G, varargin)
% ============================================================
%  RUC_FE_2D_HOMOGENEOUS.m
%
%  Callable-function refactor of RUC_FE_2D.m (delivered earlier in this
%  project). Same T6/Neo-Hookean/KUBC/Newton-Raphson-with-numerical-FD-
%  tangent FE solver, same output conventions -- but takes (H,G) directly
%  instead of a raw FALLBACK_EPS + amplitude-limit table, does not write
%  a CSV or plot a figure, and returns a results table instead.
%
%  PURPOSE IN THIS PROJECT: this is the "homogeneous, 1 RUC, kappa = 1"
%  reference solve used by compare_case_vs_homogeneous.m to normalize
%  generated NC_k####_FINAL.csv case data. With kappaMat = 1 the
%  inclusion and matrix share one Neo-Hookean material, so the body is
%  materially homogeneous -- but because the KUBC boundary field is
%  quadratic (G != 0) rather than affine, the interior field is still a
%  genuine nonlinear-equilibrium problem, NOT simply F(X) = I+H+G1*X+G2*Y
%  (see this project's own diagnosis: a naive polynomial trial field is
%  not a div(P)=0 solution once G != 0). A real FE solve is required even
%  for this "homogeneous" reference -- that is the whole reason this
%  function exists rather than a closed-form evaluation.
%
%  H, G ARE TAKEN AS GIVEN, NOT REDERIVED FROM RAW SAMPLE-SPHERE
%  COEFFICIENTS: this project's own memory notes that, under KUBC,
%  macroscopic Fbar = I + t*H holds EXACTLY (divergence theorem),
%  independent of microstructure/kappa -- so H (and, to the accuracy of
%  the Method-B boundary integration already used throughout this
%  project, G) can be read directly off any generated case's own
%  boundary F11..F22/G111..G222 columns, at any kappa. That is exactly
%  what compare_case_vs_homogeneous.m does before calling this function;
%  it means this function never needs pann_loadcases.py's sample point
%  file to reproduce a case's load direction.
%
%  USAGE
%    [T, diagOut] = RUC_FE_2D_homogeneous(H, G)
%    [T, diagOut] = RUC_FE_2D_homogeneous(H, G, 'Lx_tot',0.3, 'Ly_tot',0.3, ...)
%
%  REQUIRED
%    H : 2x2         (Fbar(t=1) = I + H)
%    G : 2x2x2       G(i,j,k), symmetric in (j,k)
%
%  NAME-VALUE OPTIONS (defaults reproduce the reference used for
%  normalization: homogeneous material, ONE repeated unit cell):
%    'Lx_tot'      (0.3)     'Ly_tot'      (0.3)
%    'NX'          (1)       'NY'          (1)
%    'zeta1'       (0.0)     'zeta2'       (0.0)
%    'Rfrac'       (0.30)    'meshFrac'    (0.03)
%    'kappaMat'    (1.0)     -- 1.0 = homogeneous (matrix==inclusion)
%    'C10_Matrix'  (400.0)   'D1_Matrix'   (0.0011429)
%    'THICKNESS'   (1.0)     'NFRAMES'     (10)
%    'relTol'      (1e-8)    'maxIter'     (30)   'maxSub' (6)
%    'h_fd'        (1e-7)    'h_energy'    (1e-6)
%    'Verbose'     (false)   'CheckAdmissibility' (true)
%    'MeshType'    ('auto')  -- 'auto' picks a pure-MATLAB regular T6
%                                grid when kappaMat==1 (no inclusion
%                                geometry, avoids the small/skewed
%                                PDE-Toolbox elements near the circle
%                                boundary that can go singular under the
%                                quadratic KUBC field); 'inclusion' forces
%                                the original PDE-Toolbox circle mesh;
%                                'regular' forces the plain grid regardless
%                                of kappaMat.
%
%  RETURNS
%    T : table, one row per frame, columns named to match the generated
%        case CSVs so ratios can be computed by simple elementwise
%        division against T_case(:, sameCols):
%          FrameID, StepTime,
%          F11,F12,F21,F22, G111,G112,G121,G122,G211,G212,G221,G222,
%          P11,P12,P21,P22, Q111,Q112,Q121,Q122,Q211,Q212,Q221,Q222,
%          P11_bnd,P12_bnd,P21_bnd,P22_bnd,
%          W_mean, ALLSE_mean  (both set from the FE's domain-integrated
%                                strain-energy density; the case CSV's
%                                two energy columns are historically
%                                distinct definitions -- see NOTE below),
%          Jmin, HillResid, QSymResid, KubcErrF, KubcErrG,
%          EnergyStressResid, ForceBalanceResid
%
%    diagOut : struct with Nnodes, Nelem, solveSeconds
%
%  NOTE ON ENERGY COLUMNS: the generated-data pipeline's ALLSE_mean and
%  W_mean come from two distinct Abaqus/Python quantities that are not
%  reproduced separately by this MATLAB solver. Both are populated here
%  from the same FE domain-integrated strain-energy density (W_density
%  in RUC_FE_2D.m's own terms) so that compare_case_vs_homogeneous.m can
%  form an energy ratio against whichever of the two the user prefers;
%  treat the resulting "energy ratio" as approximate in a way the
%  stress/gradient ratios are not.
%
%  Not executed in this environment (no MATLAB here) -- structurally
%  reviewed against the original RUC_FE_2D.m only. Run it and report any
%  error, same caveat as the original delivery.
% ============================================================

p = inputParser;
p.addParameter('Lx_tot', 0.3);
p.addParameter('Ly_tot', 0.3);
p.addParameter('NX', 1);
p.addParameter('NY', 1);
p.addParameter('zeta1', 0.0);
p.addParameter('zeta2', 0.0);
p.addParameter('Rfrac', 0.30);
p.addParameter('meshFrac', 0.05);
p.addParameter('kappaMat', 1.0);
p.addParameter('C10_Matrix', 400.0);
p.addParameter('D1_Matrix', 0.0011429);
p.addParameter('THICKNESS', 1.0);
p.addParameter('NFRAMES', 10);
p.addParameter('relTol', 1e-8);
p.addParameter('maxIter', 30);
p.addParameter('maxSub', 8);   % 8 => finest sub-step 1/256; only failing cases pay for it
p.addParameter('h_fd', 1e-7);
p.addParameter('h_energy', 1e-6);
p.addParameter('Verbose', false);
p.addParameter('CheckAdmissibility', true);
p.addParameter('MeshType', 'auto');  % 'auto' | 'regular' | 'inclusion'
p.parse(varargin{:});
o = p.Results;

tStart = tic;

% Silence the near-singular-tangent warnings that can fire transiently on
% a rejected trial Newton iterate (the backtracking line search below
% discards those iterates anyway); restored automatically on return via
% onCleanup, so this never leaks into the caller's warning state.
wState1 = warning('off', 'MATLAB:singularMatrix');
wState2 = warning('off', 'MATLAB:nearlySingularMatrix');
cleanupWarn1 = onCleanup(@() warning(wState1)); %#ok<NASGU>
cleanupWarn2 = onCleanup(@() warning(wState2)); %#ok<NASGU>

Lx_tot = o.Lx_tot; Ly_tot = o.Ly_tot;
NX = o.NX; NY = o.NY;
Lx = Lx_tot/NX; Ly = Ly_tot/NY;
Rfrac = o.Rfrac; meshFrac = o.meshFrac;
R = Rfrac*min(Lx,Ly);
meshSize = meshFrac*min(Lx,Ly);
X_CENTER = 0.5*Lx_tot; Y_CENTER = 0.5*Ly_tot;
THICKNESS = o.THICKNESS;
V0 = Lx_tot*Ly_tot*THICKNESS;
NFRAMES = o.NFRAMES;

C10_Matrix = o.C10_Matrix; D1_Matrix = o.D1_Matrix;
kappaMat = o.kappaMat;
C10_Inc = C10_Matrix*kappaMat;
D1_Inc  = D1_Matrix/kappaMat;

if o.Verbose
    fprintf('RUC_FE_2D_homogeneous: kappaMat=%.4g  Lx_tot=%.4g Ly_tot=%.4g  NX=%d NY=%d  NFRAMES=%d\n', ...
        kappaMat, Lx_tot, Ly_tot, NX, NY, NFRAMES);
end

if o.CheckAdmissibility
    checkBoundaryAdmissibility(H, G, Lx_tot, Ly_tot, X_CENTER, Y_CENTER, o.Verbose);
end

%% ---- geometry + mesh ----
switch lower(o.MeshType)
    case 'auto'
        useRegular = (kappaMat == 1);
    case 'regular'
        useRegular = true;
    case 'inclusion'
        useRegular = false;
    otherwise
        error('RUC_FE_2D_homogeneous:BadMeshType', ...
            'Unknown MeshType ''%s'' (expected auto|regular|inclusion).', o.MeshType);
end

if useRegular
    % kappaMat == 1 => matrix and "inclusion" share one material, so the
    % circle geometry is physically irrelevant here. A plain structured
    % T6 grid over the rectangle avoids the small/skewed PDE-Toolbox
    % elements near the circle-tangent points that were the first to
    % invert under the quadratic KUBC field (the reported singular-
    % tangent warning) -- and needs no toolbox.
    [nodes, elems6] = buildRegularT6Mesh(Lx_tot, Ly_tot, meshSize);
    isInclusionElem = false(size(elems6,1),1);
    if o.Verbose
        fprintf('  Mesh: regular T6 grid (kappaMat=1, no inclusion geometry)\n');
    end
else
    centers = buildInclusionCenters(NX, NY, Lx, Ly, o.zeta1, o.zeta2, Lx_tot, Ly_tot, R);
    [nodes, elems6, isInclusionElem] = buildRUCMesh(Lx_tot, Ly_tot, centers, R, meshSize);
end
Nnodes = size(nodes,1); Nelem = size(elems6,1);
if o.Verbose
    fprintf('  Mesh: %d nodes, %d T6 elements (%d inclusion, %d matrix)\n', ...
        Nnodes, Nelem, nnz(isInclusionElem), nnz(~isInclusionElem));
end

C10e = zeros(Nelem,1); D1e = zeros(Nelem,1);
C10e(isInclusionElem)  = C10_Inc;    D1e(isInclusionElem)  = D1_Inc;
C10e(~isInclusionElem) = C10_Matrix; D1e(~isInclusionElem) = D1_Matrix;

%% ---- precompute reference-config element data ----
[GP_L1,GP_L2,GP_W] = gauss3();
elemData(Nelem) = struct('dNdX',[],'detJ0',[],'w',[],'Xgp',[],'Ygp',[]);
for e = 1:Nelem
    coords6 = nodes(elems6(e,:),:);
    elemData(e).dNdX = cell(1,3);
    elemData(e).detJ0 = zeros(1,3);
    elemData(e).w = GP_W;
    elemData(e).Xgp = zeros(1,3);
    elemData(e).Ygp = zeros(1,3);
    for g = 1:3
        [~,dNdX,detJ,Xg,Yg] = t6ShapeDeriv(GP_L1(g),GP_L2(g),coords6);
        if detJ <= 0
            error('RUC_FE_2D_homogeneous:BadMesh', ...
                'Element %d has non-positive reference Jacobian at GP %d (detJ=%.4g).', e,g,detJ);
        end
        elemData(e).dNdX{g} = dNdX;
        elemData(e).detJ0(g) = detJ;
        elemData(e).Xgp(g) = Xg;
        elemData(e).Ygp(g) = Yg;
    end
end

%% ---- boundary node sets ----
TOL = 1e-6*max(Lx_tot,Ly_tot);
[leftN,rightN,botN,topN,bndAll] = classifyBoundary(nodes,Lx_tot,Ly_tot,TOL);

ndof = 2*Nnodes;
isDirichlet = false(ndof,1);
isDirichlet(2*bndAll-1) = true;
isDirichlet(2*bndAll)   = true;
freeDofs = find(~isDirichlet);

uBC1 = buildDirichletTarget(nodes, isDirichlet, H,G, X_CENTER,Y_CENTER);

%% ---- load stepping / Newton-Raphson ----
u = zeros(ndof,1);
tNow = 0.0;
tlist = linspace(1/NFRAMES,1,NFRAMES);
results = struct([]);

for iFrame = 1:NFRAMES
    tTarget = tlist(iFrame);
    [u, ok] = advanceToTarget(u, tNow, tTarget, elems6, elemData, C10e, D1e, ...
        isDirichlet, freeDofs, uBC1, THICKNESS, o.h_fd, o.relTol, o.maxIter, o.maxSub);
    if ~ok
        error('RUC_FE_2D_homogeneous:NoConverge', ...
            'Frame %d (t=%.4f) failed to converge even after adaptive sub-stepping.', iFrame, tTarget);
    end
    tNow = tTarget;
    res = postprocessFrame(u, nodes, elems6, elemData, C10e, D1e, ...
        leftN,rightN,botN,topN,bndAll, X_CENTER,Y_CENTER, Lx_tot,Ly_tot, THICKNESS, V0, H,G,tTarget, o.h_energy);
    res.FrameID = iFrame;
    res.StepTime = tTarget;
    if res.Jmin <= 0
        error('RUC_FE_2D_homogeneous:InvertedElement', ...
            ['Frame %d (t=%.4f) satisfied the Newton residual tolerance but contains ' ...
             'an inverted/degenerate element (Jmin=%.4g <= 0) -- not a physically ' ...
             'admissible solution. Treating this as a failed frame.'], iFrame, tTarget, res.Jmin);
    end
    if isempty(results), results = res; else, results(end+1) = res; end %#ok<AGROW>
    if o.Verbose
        fprintf('  Frame %2d  t=%.3f  Jmin=%7.4f  HillResid=%.2e  QSymResid=%.2e  KubcErrF=%.2e  KubcErrG=%.2e\n', ...
            iFrame, tTarget, res.Jmin, res.HillResid, res.QSymResid, res.KubcErrF, res.KubcErrG);
    end
end

%% ---- assemble output table with case-CSV-compatible column names ----
Traw = struct2table(results);
T = table();
T.FrameID  = Traw.FrameID;
T.StepTime = Traw.StepTime;
T.F11 = Traw.F11_bnd; T.F12 = Traw.F12_bnd; T.F21 = Traw.F21_bnd; T.F22 = Traw.F22_bnd;
T.G111 = Traw.G111; T.G112 = Traw.G112; T.G121 = Traw.G121; T.G122 = Traw.G122;
T.G211 = Traw.G211; T.G212 = Traw.G212; T.G221 = Traw.G221; T.G222 = Traw.G222;
T.P11 = Traw.P11; T.P12 = Traw.P12; T.P21 = Traw.P21; T.P22 = Traw.P22;
T.Q111 = Traw.Q111; T.Q112 = Traw.Q112; T.Q121 = Traw.Q121; T.Q122 = Traw.Q122;
T.Q211 = Traw.Q211; T.Q212 = Traw.Q212; T.Q221 = Traw.Q221; T.Q222 = Traw.Q222;
T.P11_bnd = Traw.P11_bnd; T.P12_bnd = Traw.P12_bnd; T.P21_bnd = Traw.P21_bnd; T.P22_bnd = Traw.P22_bnd;
T.W_mean     = Traw.W_density;   % see NOTE in header
T.ALLSE_mean = Traw.W_density;   % see NOTE in header
T.Jmin = Traw.Jmin;
T.HillResid = Traw.HillResid;
T.QSymResid = Traw.QSymResid;
T.KubcErrF = Traw.KubcErrF;
T.KubcErrG = Traw.KubcErrG;
T.EnergyStressResid = Traw.EnergyStressResid;
T.ForceBalanceResid = Traw.ForceBalanceResid;

diagOut = struct('Nnodes',Nnodes,'Nelem',Nelem,'solveSeconds',toc(tStart));
if o.Verbose
    fprintf('  done in %.1f s\n', diagOut.solveSeconds);
end

end % RUC_FE_2D_homogeneous


% ============================================================
% LOCAL FUNCTIONS -- unchanged from RUC_FE_2D.m (see that file for the
% full derivations/comments); only checkBoundaryAdmissibility gained a
% 'verbose' argument so it can stay silent inside a 273-case loop.
% ============================================================

function F = FfromPoly(H,G,X,Y,Xc,Yc)
    xr = X-Xc; yr = Y-Yc;
    F = eye(2);
    for i = 1:2
        for j = 1:2
            F(i,j) = F(i,j) + H(i,j) + G(i,j,1)*xr + G(i,j,2)*yr;
        end
    end
end

function checkBoundaryAdmissibility(H,G,Lx_tot,Ly_tot,Xc,Yc,verbose)
    ns = 200;
    Xs = linspace(0,Lx_tot,ns); Ys = linspace(0,Ly_tot,ns);
    worstDet = inf;
    for x = Xs
        worstDet = min(worstDet, det(FfromPoly(H,G,x,0,Xc,Yc)));
        worstDet = min(worstDet, det(FfromPoly(H,G,x,Ly_tot,Xc,Yc)));
    end
    for y = Ys
        worstDet = min(worstDet, det(FfromPoly(H,G,0,y,Xc,Yc)));
        worstDet = min(worstDet, det(FfromPoly(H,G,Lx_tot,y,Xc,Yc)));
    end
    if worstDet <= 0
        warning(['KUBC boundary polynomial is INADMISSIBLE at t=1: min det(F) on ' ...
                 'the boundary = %.4g <= 0.'], worstDet);
    elseif verbose
        fprintf('  Boundary admissibility OK: min det(F) on boundary at t=1 = %.4g > 0\n', worstDet);
    end
end

function u = boundaryDispVal(H,G,X,Y,Xc,Yc,t)
    xr = X-Xc; yr = Y-Yc;
    u = zeros(2,1);
    for i = 1:2
        ui = H(i,1)*xr + H(i,2)*yr ...
           + 0.5*G(i,1,1)*xr^2 + G(i,1,2)*xr*yr + 0.5*G(i,2,2)*yr^2;
        u(i) = t*ui;
    end
end

function uBC1 = buildDirichletTarget(nodes, isDirichlet, H,G,Xc,Yc)
    ndof = 2*size(nodes,1);
    uBC1 = zeros(ndof,1);
    for n = 1:size(nodes,1)
        dx = 2*n-1; dy = 2*n;
        if isDirichlet(dx) || isDirichlet(dy)
            uval = boundaryDispVal(H,G,nodes(n,1),nodes(n,2),Xc,Yc,1.0);
            uBC1(dx) = uval(1);
            uBC1(dy) = uval(2);
        end
    end
end

function tf = circleIntersectsDomain(cx,cy,Rc,xmin,xmax,ymin,ymax)
    tf = ~(cx+Rc < xmin || cx-Rc > xmax || cy+Rc < ymin || cy-Rc > ymax);
end

function centers = buildInclusionCenters(NX,NY,Lx,Ly,zeta1,zeta2,Lx_tot,Ly_tot,R)
    centers = zeros(0,2);
    for i = -1:NX
        for j = -1:NY
            cx = (i+0.5)*Lx + zeta1;
            cy = (j+0.5)*Ly + zeta2;
            if circleIntersectsDomain(cx,cy,R,0,Lx_tot,0,Ly_tot)
                centers(end+1,:) = [cx,cy]; %#ok<AGROW>
            end
        end
    end
    if isempty(centers)
        error('RUC_FE_2D_homogeneous:NoInclusions', ...
            'No inclusion circles intersect the domain -- check NX/NY/zeta/geometry.');
    end
end

function [nodes, elems6, isInclusionElem] = buildRUCMesh(Lx_tot,Ly_tot,centers,R,meshSize)
    if exist('createpde','file') ~= 2
        error('RUC_FE_2D_homogeneous:NoPDEToolbox', ['This mesh generator needs the Partial ' ...
            'Differential Equation Toolbox (createpde/geometryFromEdges/generateMesh).']);
    end
    Nc = size(centers,1);
    rectGd = [3;4; 0;Lx_tot;Lx_tot;0; 0;0;Ly_tot;Ly_tot];
    gd = rectGd;
    names = {'R1'};
    sf = 'R1';
    for k = 1:Nc
        circGd = [1; centers(k,1); centers(k,2); R; 0;0;0;0;0;0];
        gd = [gd, circGd]; %#ok<AGROW>
        nm = sprintf('C%d',k);
        names{end+1} = nm; %#ok<AGROW>
        sf = [sf '+' nm]; %#ok<AGROW>
    end
    ns = char(names)';
    dl = decsg(gd, sf, ns);

    model = createpde(1);
    geometryFromEdges(model, dl);
    generateMesh(model, 'Hmax', meshSize, 'GeometricOrder', 'quadratic');

    nodes = model.Mesh.Nodes';
    elems6 = model.Mesh.Elements';

    Nelem = size(elems6,1);
    isInclusionElem = false(Nelem,1);
    for e = 1:Nelem
        cn = elems6(e,1:3);
        cx = mean(nodes(cn,1)); cy = mean(nodes(cn,2));
        d2 = (centers(:,1)-cx).^2 + (centers(:,2)-cy).^2;
        isInclusionElem(e) = any(d2 <= R^2);
    end
end

function [nodes, elems6] = buildRegularT6Mesh(Lx_tot, Ly_tot, meshSize)
% Pure-MATLAB structured grid of CCW right-triangle T6 elements over the
% plain rectangle [0,Lx_tot]x[0,Ly_tot]. No PDE Toolbox dependency. Used
% (via 'MeshType','auto') whenever kappaMat==1, where the inclusion
% geometry is physically irrelevant and the PDE-Toolbox circle-refined
% mesh's small/skewed elements near the circle-tangent points were the
% first to invert under the quadratic KUBC field.
    nX = max(1, round(Lx_tot/meshSize));
    nY = max(1, round(Ly_tot/meshSize));
    dx = Lx_tot/nX; dy = Ly_tot/nY;

    [Xc, Yc] = meshgrid((0:nX)*dx, (0:nY)*dy);
    cornerID = reshape(1:numel(Xc), nY+1, nX+1);
    nCorner = numel(Xc);

    [Xh, Yh] = meshgrid(((0:nX-1)+0.5)*dx, (0:nY)*dy);
    hmidID = nCorner + reshape(1:numel(Xh), nY+1, nX);
    nHmid = numel(Xh);

    [Xv, Yv] = meshgrid((0:nX)*dx, ((0:nY-1)+0.5)*dy);
    vmidID = nCorner + nHmid + reshape(1:numel(Xv), nY, nX+1);
    nVmid = numel(Xv);

    [Xd, Yd] = meshgrid(((0:nX-1)+0.5)*dx, ((0:nY-1)+0.5)*dy);
    diagID = nCorner + nHmid + nVmid + reshape(1:numel(Xd), nY, nX);

    nodes = [Xc(:) Yc(:); Xh(:) Yh(:); Xv(:) Yv(:); Xd(:) Yd(:)];

    elems6 = zeros(2*nX*nY, 6);
    e = 0;
    for i = 0:nX-1
        for j = 0:nY-1
            c00 = cornerID(j+1, i+1);
            c10 = cornerID(j+1, i+2);
            c11 = cornerID(j+2, i+2);
            c01 = cornerID(j+2, i+1);
            hBot = hmidID(j+1, i+1);
            hTop = hmidID(j+2, i+1);
            vLef = vmidID(j+1, i+1);
            vRig = vmidID(j+1, i+2);
            dMid = diagID(j+1, i+1);

            % Two CCW right triangles per rectangular cell, split along
            % the c00-c11 diagonal, each with its own edge-midpoint node
            % (vertex order matches t6ShapeDeriv's L1,L2,L3,mid12,mid23,
            % mid31 convention).
            e = e+1;
            elems6(e,:) = [c00, c10, c11, hBot, vRig, dMid];
            e = e+1;
            elems6(e,:) = [c00, c11, c01, dMid, hTop, vLef];
        end
    end
end

function [leftN,rightN,botN,topN,bndAll] = classifyBoundary(nodes,Lx,Ly,TOL)
    x = nodes(:,1); y = nodes(:,2);
    leftN  = find(abs(x-0)  < TOL);
    rightN = find(abs(x-Lx) < TOL);
    botN   = find(abs(y-0)  < TOL);
    topN   = find(abs(y-Ly) < TOL);
    bndAll = unique([leftN;rightN;botN;topN]);
end

function [L1w,L2w,W] = gauss3()
    L1w = [2/3, 1/6, 1/6];
    L2w = [1/6, 2/3, 1/6];
    W   = [1/6, 1/6, 1/6];
end

function [N, dNdX, detJ, X, Y] = t6ShapeDeriv(L1,L2,coords6)
    L3 = 1-L1-L2;
    N = [L1*(2*L1-1); L2*(2*L2-1); L3*(2*L3-1); 4*L1*L2; 4*L2*L3; 4*L3*L1];
    dN_dL1 = [4*L1-1; 0; -(4*L3-1); 4*L2; -4*L2; 4*(L3-L1)];
    dN_dL2 = [0; 4*L2-1; -(4*L3-1); 4*L1; 4*(L3-L2); -4*L1];
    x = coords6(:,1); y = coords6(:,2);
    dXdL1 = dN_dL1'*x; dXdL2 = dN_dL2'*x;
    dYdL1 = dN_dL1'*y; dYdL2 = dN_dL2'*y;
    J = [dXdL1 dXdL2; dYdL1 dYdL2];
    detJ = det(J);
    Jinv = inv(J);
    dL1dX = Jinv(1,1); dL2dX = Jinv(2,1);
    dL1dY = Jinv(1,2); dL2dY = Jinv(2,2);
    dNdX = zeros(6,2);
    dNdX(:,1) = dN_dL1*dL1dX + dN_dL2*dL2dX;
    dNdX(:,2) = dN_dL1*dL1dY + dN_dL2*dL2dY;
    X = N'*x; Y = N'*y;
end

function W = neoHookeanW(F, C10, D1)
    J = det(F);
    I1 = F(1,1)^2+F(2,1)^2 + F(1,2)^2+F(2,2)^2 + 1;
    W = C10*(J^(-2/3)*I1 - 3) + (1/D1)*(J-1)^2;
end

function P = neoHookeanP(F, C10, D1)
    J = det(F);
    Finv = inv(F);
    FinvT = Finv';
    I1 = F(1,1)^2+F(2,1)^2 + F(1,2)^2+F(2,2)^2 + 1;
    P = C10*J^(-2/3)*(2*F - (2/3)*I1*FinvT) + (2/D1)*(J-1)*J*FinvT;
end

function Pfd = neoHookeanP_FD(F, C10, D1, h)
    Pfd = zeros(2,2);
    for i = 1:2
        for j = 1:2
            Fp = F; Fp(i,j) = Fp(i,j) + h;
            Fm = F; Fm(i,j) = Fm(i,j) - h;
            Wp = neoHookeanW(Fp, C10, D1);
            Wm = neoHookeanW(Fm, C10, D1);
            Pfd(i,j) = (Wp - Wm) / (2*h);
        end
    end
end

function Re = elementResidual(ue, ed, C10, D1, THICKNESS)
    u1 = ue(1:2:end); u2 = ue(2:2:end);
    Re = zeros(12,1);
    for g = 1:3
        dNdX = ed.dNdX{g};
        F = eye(2);
        F(1,1) = F(1,1) + dNdX(:,1)'*u1;  F(1,2) = F(1,2) + dNdX(:,2)'*u1;
        F(2,1) = F(2,1) + dNdX(:,1)'*u2;  F(2,2) = F(2,2) + dNdX(:,2)'*u2;
        P = neoHookeanP(F, C10, D1);
        wdv = ed.w(g)*ed.detJ0(g)*THICKNESS;
        for a = 1:6
            fa = P*dNdX(a,:)';
            Re(2*a-1) = Re(2*a-1) + wdv*fa(1);
            Re(2*a)   = Re(2*a)   + wdv*fa(2);
        end
    end
end

function Ke = elementTangentFD(ue, ed, C10, D1, THICKNESS, h)
    Ke = zeros(12,12);
    for j = 1:12
        up = ue; up(j) = up(j)+h;
        um = ue; um(j) = um(j)-h;
        Rp = elementResidual(up, ed, C10, D1, THICKNESS);
        Rm = elementResidual(um, ed, C10, D1, THICKNESS);
        Ke(:,j) = (Rp-Rm)/(2*h);
    end
end

function [Rfull, Kff] = assembleResidualTangent(u, elems6, elemData, C10e, D1e, THICKNESS, freeDofs, h_fd)
    Nelem = size(elems6,1);
    ndof = numel(u);
    Rfull = zeros(ndof,1);
    Ii = zeros(Nelem*144,1); Jj = zeros(Nelem*144,1); Vv = zeros(Nelem*144,1);
    ptr = 0;
    for e = 1:Nelem
        nodesE = elems6(e,:);
        gdofs = zeros(12,1);
        gdofs(1:2:end) = 2*nodesE-1; gdofs(2:2:end) = 2*nodesE;
        ue = u(gdofs);
        Re = elementResidual(ue, elemData(e), C10e(e), D1e(e), THICKNESS);
        Rfull(gdofs) = Rfull(gdofs) + Re;
        Ke = elementTangentFD(ue, elemData(e), C10e(e), D1e(e), THICKNESS, h_fd);
        for a = 1:12
            for b = 1:12
                ptr = ptr+1;
                Ii(ptr) = gdofs(a); Jj(ptr) = gdofs(b); Vv(ptr) = Ke(a,b);
            end
        end
    end
    Kfull = sparse(Ii(1:ptr), Jj(1:ptr), Vv(1:ptr), ndof, ndof);
    Kff = Kfull(freeDofs, freeDofs);
end

function Rfull = assembleResidualOnly(u, elems6, elemData, C10e, D1e, THICKNESS)
    Nelem = size(elems6,1);
    Rfull = zeros(numel(u),1);
    for e = 1:Nelem
        nodesE = elems6(e,:);
        gdofs = zeros(12,1);
        gdofs(1:2:end) = 2*nodesE-1; gdofs(2:2:end) = 2*nodesE;
        ue = u(gdofs);
        Re = elementResidual(ue, elemData(e), C10e(e), D1e(e), THICKNESS);
        Rfull(gdofs) = Rfull(gdofs) + Re;
    end
end

function [u, converged, resNorm] = solveLoadStep(u0, tTarget, elems6, elemData, C10e, D1e, ...
        isDirichlet, freeDofs, uBC1, THICKNESS, h_fd, relTol, maxIter)
    u = u0;
    u(isDirichlet) = tTarget*uBC1(isDirichlet);
    converged = false;
    R0norm = [];
    resNorm = inf;
    for it = 1:maxIter
        [Rfull, Kff] = assembleResidualTangent(u, elems6, elemData, C10e, D1e, THICKNESS, freeDofs, h_fd);
        resNorm = norm(Rfull(freeDofs));
        if isempty(R0norm), R0norm = max(resNorm, 1e-30); end
        if resNorm < max(1e-12, relTol*R0norm)
            converged = true; break;
        end
        du = -(Kff \ Rfull(freeDofs));

        % Backtracking line search: a full Newton step can walk a trial
        % iterate through a near-singular/inverted deformation gradient
        % (det(F)<=0), which is exactly what produced the reported
        % "Matrix is singular" warning out of neoHookeanP's inv(F). Halve
        % the step until the trial residual is finite and no worse than
        % the current one, rather than accepting an ill-conditioned step
        % outright.
        alpha = 1.0;
        uTrial = u; uTrial(freeDofs) = u(freeDofs) + alpha*du;
        RTrial = assembleResidualOnly(uTrial, elems6, elemData, C10e, D1e, THICKNESS);
        resTrialNorm = norm(RTrial(freeDofs));
        nBacktrack = 0;
        while (~isfinite(resTrialNorm) || resTrialNorm > 1.0001*resNorm) && nBacktrack < 8
            alpha = 0.5*alpha;
            uTrial = u; uTrial(freeDofs) = u(freeDofs) + alpha*du;
            RTrial = assembleResidualOnly(uTrial, elems6, elemData, C10e, D1e, THICKNESS);
            resTrialNorm = norm(RTrial(freeDofs));
            nBacktrack = nBacktrack + 1;
        end
        u = uTrial;
    end
end

function [u, converged] = advanceToTarget(u0, tCur, tTarget, elems6, elemData, C10e, D1e, ...
        isDirichlet, freeDofs, uBC1, THICKNESS, h_fd, relTol, maxIter, maxSub)
    u = u0; tNow = tCur;
    subTargets = tTarget;
    nSub = 0;
    while ~isempty(subTargets)
        tTry = subTargets(1);
        [uTry, conv] = solveLoadStep(u, tTry, elems6, elemData, C10e, D1e, ...
            isDirichlet, freeDofs, uBC1, THICKNESS, h_fd, relTol, maxIter);
        if conv
            u = uTry; tNow = tTry; %#ok<NASGU>
            subTargets(1) = [];
        else
            nSub = nSub+1;
            if nSub > maxSub
                converged = false; return;
            end
            tMid = 0.5*(tNow+tTry);
            subTargets = [tMid, subTargets]; %#ok<AGROW>
        end
    end
    converged = true;
end

function w = trapWeights1D(p)
    n = numel(p);
    w = zeros(n,1);
    if n == 1, return; end
    w(1) = 0.5*(p(2)-p(1));
    w(n) = 0.5*(p(n)-p(n-1));
    for i = 2:n-1
        w(i) = 0.5*(p(i+1)-p(i-1));
    end
end

function [dat,w] = edgeData(u,nodes,nlist,sortCol)
    x = nodes(nlist,1); y = nodes(nlist,2);
    u1 = u(2*nlist-1); u2 = u(2*nlist);
    dat = [x,y,u1,u2];
    [~,ord] = sort(dat(:,sortCol));
    dat = dat(ord,:);
    w = trapWeights1D(dat(:,sortCol));
end

function [Fbar,Gbar] = boundaryFbarG_FE_matlab(u, nodes, leftN,rightN,botN,topN, Lx_tot,Ly_tot, Xc,Yc)
    Vtot = Lx_tot*Ly_tot;
    [left,wL]   = edgeData(u,nodes,leftN, 2);
    [right,wR]  = edgeData(u,nodes,rightN,2);
    [bottom,wB] = edgeData(u,nodes,botN,  1);
    [top,wT]    = edgeData(u,nodes,topN,  1);

    u1L=left(:,3); u2L=left(:,4); ycL=left(:,2)-Yc;
    u1R=right(:,3);u2R=right(:,4);ycR=right(:,2)-Yc;
    u1B=bottom(:,3);u2B=bottom(:,4);xcB=bottom(:,1)-Xc;
    u1T=top(:,3);   u2T=top(:,4);   xcT=top(:,1)-Xc;

    L_Iu1=wL'*u1L; L_Iu2=wL'*u2L;
    R_Iu1=wR'*u1R; R_Iu2=wR'*u2R;
    B_Iu1=wB'*u1B; B_Iu2=wB'*u2B;
    T_Iu1=wT'*u1T; T_Iu2=wT'*u2T;

    xcL_c = left(1,1)-Xc;  xcR_c = right(1,1)-Xc;
    ycB_c = bottom(1,2)-Yc; ycT_c = top(1,2)-Yc;

    L_Iu1x1=xcL_c*L_Iu1; L_Iu1x2=wL'*(u1L.*ycL);
    L_Iu2x1=xcL_c*L_Iu2; L_Iu2x2=wL'*(u2L.*ycL);
    R_Iu1x1=xcR_c*R_Iu1; R_Iu1x2=wR'*(u1R.*ycR);
    R_Iu2x1=xcR_c*R_Iu2; R_Iu2x2=wR'*(u2R.*ycR);
    B_Iu1x1=wB'*(u1B.*xcB); B_Iu1x2=ycB_c*B_Iu1;
    B_Iu2x1=wB'*(u2B.*xcB); B_Iu2x2=ycB_c*B_Iu2;
    T_Iu1x1=wT'*(u1T.*xcT); T_Iu1x2=ycT_c*T_Iu1;
    T_Iu2x1=wT'*(u2T.*xcT); T_Iu2x2=ycT_c*T_Iu2;

    F11 = 1 + (-L_Iu1+R_Iu1)/Vtot;
    F12 = (-B_Iu1+T_Iu1)/Vtot;
    F21 = (-L_Iu2+R_Iu2)/Vtot;
    F22 = 1 + (-B_Iu2+T_Iu2)/Vtot;

    S111=(-L_Iu1x1+R_Iu1x1)/Vtot; S112=(-L_Iu1x2+R_Iu1x2)/Vtot;
    S121=(-B_Iu1x1+T_Iu1x1)/Vtot; S122=(-B_Iu1x2+T_Iu1x2)/Vtot;
    S211=(-L_Iu2x1+R_Iu2x1)/Vtot; S212=(-L_Iu2x2+R_Iu2x2)/Vtot;
    S221=(-B_Iu2x1+T_Iu2x1)/Vtot; S222=(-B_Iu2x2+T_Iu2x2)/Vtot;

    varX = Lx_tot^2/12; varY = Ly_tot^2/12;
    A1 = 0.25*(S111+S122); A2 = 0.25*(S211+S222);

    G111=(S111-A1)/varX; G112=S112/varY;
    G121=S121/varX;      G122=(S122-A1)/varY;
    G211=(S211-A2)/varX; G212=S212/varY;
    G221=S221/varX;      G222=(S222-A2)/varY;

    Fbar = [F11 F12; F21 F22];
    Gbar = zeros(2,2,2);
    Gbar(1,1,1)=G111; Gbar(1,1,2)=G112; Gbar(1,2,1)=G121; Gbar(1,2,2)=G122;
    Gbar(2,1,1)=G211; Gbar(2,1,2)=G212; Gbar(2,2,1)=G221; Gbar(2,2,2)=G222;
end

function [Fx,Fy] = sideReactionForces(Rfull, nodeList)
    Fx = sum(Rfull(2*nodeList-1));
    Fy = sum(Rfull(2*nodeList));
end

function res = postprocessFrame(u, nodes, elems6, elemData, C10e,D1e, ...
        leftN,rightN,botN,topN,bndAll, Xc,Yc, Lx_tot,Ly_tot,THICKNESS,V0, H,G,tTarget, h_energy)

    Nelem = size(elems6,1);
    Pnum = zeros(2,2);
    Qnum = zeros(2,2,2);
    Jmin = inf;
    Wtot = 0;
    EnergyStressResid = 0;
    for e = 1:Nelem
        nodesE = elems6(e,:);
        gdofs = zeros(12,1);
        gdofs(1:2:end) = 2*nodesE-1; gdofs(2:2:end) = 2*nodesE;
        ue = u(gdofs);
        u1 = ue(1:2:end); u2 = ue(2:2:end);
        for g = 1:3
            dNdX = elemData(e).dNdX{g};
            F = eye(2);
            F(1,1) = F(1,1) + dNdX(:,1)'*u1;  F(1,2) = F(1,2) + dNdX(:,2)'*u1;
            F(2,1) = F(2,1) + dNdX(:,1)'*u2;  F(2,2) = F(2,2) + dNdX(:,2)'*u2;
            Jg = det(F); Jmin = min(Jmin,Jg);
            P = neoHookeanP(F, C10e(e), D1e(e));
            wdv = elemData(e).w(g)*elemData(e).detJ0(g)*THICKNESS;
            Xr = elemData(e).Xgp(g)-Xc;
            Yr = elemData(e).Ygp(g)-Yc;
            Pnum = Pnum + wdv*P;
            Qnum(:,:,1) = Qnum(:,:,1) + wdv*Xr*P;
            Qnum(:,:,2) = Qnum(:,:,2) + wdv*Yr*P;

            Wgp = neoHookeanW(F, C10e(e), D1e(e));
            Wtot = Wtot + wdv*Wgp;
            Pfd = neoHookeanP_FD(F, C10e(e), D1e(e), h_energy);
            gpScale = max(1e-30, max(abs(P(:))));
            gpResid = max(abs(P(:)-Pfd(:))) / gpScale;
            EnergyStressResid = max(EnergyStressResid, gpResid);
        end
    end
    Pbar = Pnum/V0;
    Qbar = Qnum/V0;
    W_density = Wtot/V0;

    Rfull = assembleResidualOnly(u, elems6, elemData, C10e, D1e, THICKNESS);

    [RF_L_x, RF_L_y] = sideReactionForces(Rfull, leftN);
    [RF_R_x, RF_R_y] = sideReactionForces(Rfull, rightN);
    [RF_B_x, RF_B_y] = sideReactionForces(Rfull, botN);
    [RF_T_x, RF_T_y] = sideReactionForces(Rfull, topN);

    HillPV0 = zeros(2,2);
    QB_xxx=0; QB_xyy=0; QB_xxy=0;
    QB_yxx=0; QB_yyy=0; QB_yxy=0;
    sumRFx = 0; sumRFy = 0;
    for idx = 1:numel(bndAll)
        nlab = bndAll(idx);
        xr = nodes(nlab,1)-Xc; yr = nodes(nlab,2)-Yc;
        Fx = Rfull(2*nlab-1); Fy = Rfull(2*nlab);
        HillPV0(1,1)=HillPV0(1,1)+Fx*xr; HillPV0(1,2)=HillPV0(1,2)+Fx*yr;
        HillPV0(2,1)=HillPV0(2,1)+Fy*xr; HillPV0(2,2)=HillPV0(2,2)+Fy*yr;
        QB_xxx=QB_xxx+Fx*xr*xr; QB_xyy=QB_xyy+Fx*yr*yr; QB_xxy=QB_xxy+Fx*xr*yr;
        QB_yxx=QB_yxx+Fy*xr*xr; QB_yyy=QB_yyy+Fy*yr*yr; QB_yxy=QB_yxy+Fy*xr*yr;
        sumRFx = sumRFx + Fx; sumRFy = sumRFy + Fy;
    end
    Pbar_bnd = HillPV0/V0;
    hillScale = max(1e-30, max(abs(Pnum(:))));
    HillResid = max(abs(Pnum(:)-HillPV0(:))) / hillScale;

    forceScale = max(1e-30, max(abs([RF_L_x,RF_L_y,RF_R_x,RF_R_y,RF_B_x,RF_B_y,RF_T_x,RF_T_y])));
    ForceBalanceResid = norm([sumRFx,sumRFy]) / forceScale;

    domSym_1_11=2*Qnum(1,1,1); domSym_1_22=2*Qnum(1,2,2); domSym_1_12=Qnum(1,1,2)+Qnum(1,2,1);
    domSym_2_11=2*Qnum(2,1,1); domSym_2_22=2*Qnum(2,2,2); domSym_2_12=Qnum(2,1,2)+Qnum(2,2,1);
    qDiffs = [domSym_1_11-QB_xxx, domSym_1_22-QB_xyy, domSym_1_12-QB_xxy, ...
              domSym_2_11-QB_yxx, domSym_2_22-QB_yyy, domSym_2_12-QB_yxy];
    qScale = max(1e-30, max(abs(Qnum(:))));
    QSymResid = max(abs(qDiffs)) / qScale;

    [Fbar_bnd, Gbar_bnd] = boundaryFbarG_FE_matlab(u, nodes, leftN,rightN,botN,topN, Lx_tot,Ly_tot, Xc,Yc);

    KubcErrF = max(max(abs(Fbar_bnd - (eye(2)+tTarget*H))));
    KubcErrG = 0;
    for i = 1:2
        for j = 1:2
            for k = 1:2
                KubcErrG = max(KubcErrG, abs(Gbar_bnd(i,j,k)-tTarget*G(i,j,k)));
            end
        end
    end

    res = struct();
    res.F11_bnd=Fbar_bnd(1,1); res.F12_bnd=Fbar_bnd(1,2);
    res.F21_bnd=Fbar_bnd(2,1); res.F22_bnd=Fbar_bnd(2,2);
    res.G111=Gbar_bnd(1,1,1); res.G112=Gbar_bnd(1,1,2); res.G121=Gbar_bnd(1,2,1); res.G122=Gbar_bnd(1,2,2);
    res.G211=Gbar_bnd(2,1,1); res.G212=Gbar_bnd(2,1,2); res.G221=Gbar_bnd(2,2,1); res.G222=Gbar_bnd(2,2,2);
    res.P11=Pbar(1,1); res.P12=Pbar(1,2); res.P21=Pbar(2,1); res.P22=Pbar(2,2);
    res.Q111=Qbar(1,1,1); res.Q112=Qbar(1,1,2); res.Q121=Qbar(1,2,1); res.Q122=Qbar(1,2,2);
    res.Q211=Qbar(2,1,1); res.Q212=Qbar(2,1,2); res.Q221=Qbar(2,2,1); res.Q222=Qbar(2,2,2);
    res.P11_bnd=Pbar_bnd(1,1); res.P12_bnd=Pbar_bnd(1,2);
    res.P21_bnd=Pbar_bnd(2,1); res.P22_bnd=Pbar_bnd(2,2);
    res.W_density=W_density;
    res.HillResid=HillResid; res.QSymResid=QSymResid; res.Jmin=Jmin;
    res.KubcErrF=KubcErrF; res.KubcErrG=KubcErrG;
    res.EnergyStressResid=EnergyStressResid;
    res.ForceBalanceResid=ForceBalanceResid;
end