//
// This file is part of MorphoDynamX - http://www.MorphoDynamX.org
// Copyright (C) 2012-2017 Richard S. Smith and collaborators.
//
// If you use MorphoDynamX in your work, please cite:
//   http://dx.doi.org/10.7554/eLife.05864
//
// MorphoDynamX is free software, and is licensed under under the terms of the 
// GNU General (GPL) Public License version 2.0, http://www.gnu.org/licenses.
//
#include <CellDisk.hpp>
#include <CCUtils.hpp>

namespace CellDisk
{
  /*
   * Auxin model processes
   */	
  bool AuxinGradient::initialize(CellTissue &_tissue, CCIndexDataAttr &_indexAttr, CellDataAttr &_cellAttr, EdgeDataAttr &_edgeAttr)
  {
    tissue = &_tissue;
    indexAttr = &_indexAttr;
    cellAttr = &_cellAttr;
    edgeAttr = &_edgeAttr;

    if(!tissue)
      throw(QString("AuxinGradient::initialize Tissue empty"));
    if(!indexAttr or !cellAttr or !edgeAttr)
      throw(QString("AuxinGradient::initialize Attribute map empty"));

    Prod = parm("Production").toDouble();
    Trans = parm("Transport").toDouble();
    Diff = parm("Diffusion").toDouble();
    Decay = parm("Decay").toDouble();

    return true;
  }

  void AuxinGradient::initDerivs(SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr)
  {
    tissue->updateGeometry();
    forall(CCIndex c, tissue->dualGraph().vertices()) {
      CellData &cCD = (*cellAttr)[c];
      CCIndexData &cIdx = (*indexAttr)[c];
      cCD.area = cIdx.measure;
      cCD.perimeter = 0;

      forall(const Flip &flip, neighborVertexFlips(tissue->dualGraph(), c)) {
        CCIndex w = flip.interior;

        CCIndexData &wIdx = (*indexAttr)[w];
        cCD.perimeter += wIdx.measure;
      }
    }
  }

  // Calculate derivative for the solver
  void AuxinGradient::calcDerivatives(const SolverT &solver, VertexAttr &vertexAttr)
  {
    if(!tissue) {
      printf("CellGrid::calcDerivatives Tissue is null\n");
      return;
    }
    CCStructure &cs = tissue->cellStructure();
    CCStructure &dg = tissue->dualGraph();

    // First clear derivatives
    for(CCIndex c : dg.vertices())
      vertexAttr[c].dx.zero();

    // Clear PIN
    for(auto &pr : *edgeAttr)
      pr.second.pin = 0;

    // Now calculate
    for(CCIndex c : dg.vertices()) {
      CellData &cCD = (*cellAttr)[c];
      if(cCD.fixedConc)
        cCD.a = cCD.amount;
    }
    // Now calculate
    for(CCIndex c : dg.vertices()) {
      CellData &cCD = (*cellAttr)[c];
      CCIndexData &cIdx = (*indexAttr)[c];
      SolverT::SolverVertexData &cSD = vertexAttr[c];

      // Calculate pin
      double aTotal = 0;
      for(const Flip &flip : neighborVertexFlips(tissue->dualGraph(), c)) {
        CCIndex w = flip.interior;
        CCIndex n = flip.otherFacet(c);
  
        CellData &nCD = (*cellAttr)[n];
        CCIndexData &wIdx = (*indexAttr)[w];
        aTotal += nCD.a * wIdx.measure;
      }
      // Production
      if(cCD.cellType == SourceCell)
        cSD.dx[0] += cCD.amount;
      else
        cSD.dx[0] += Prod;
  
      // Decay
      if(cCD.cellType == SinkCell)
        cSD.dx[0] -= cCD.amount * cCD.a;
      else
        cSD.dx[0] -= Decay * cCD.a;
  
      // Calculate diffusion and transport
      for(const Flip &flip : neighborVertexFlips(dg, c)) {
        CCIndex w = flip.interior;
        CCIndex n = flip.otherFacet(c);
  
        // Diffusion
        CellData &nCD = (*cellAttr)[n];
        SolverT::SolverVertexData &nSD = vertexAttr[n];
        CCIndexData &wIdx = (*indexAttr)[w];
        if (cCD.fixedConc == true){
          cSD.dx[0] = 0.;
        }
        else 
         cSD.dx[0] += Diff * (nCD.a - cCD.a) * wIdx.measure / cIdx.measure;
  
        // Transport
        if(aTotal <= 0)
          continue;
        CCSignedIndex e = orientedMembrane(cs, c, n);
        EdgeData &eED = (*edgeAttr)[e];
        eED.pin = (nCD.a * wIdx.measure)/aTotal;
        CCIndexData &nIdx = (*indexAttr)[n];
        cSD.dx[0] -= Trans * cCD.a * eED.pin / cIdx.measure;
        nSD.dx[0] += Trans * cCD.a * eED.pin / nIdx.measure;
      }
    }
  }

  // Initialize the main solver process
  bool AuxinSolver::initialize(QWidget *parent)
  {
    mesh = currentMesh();
    if(!mesh)
      throw(QString("AuxinSolver::initialize No current mesh"));
    indexAttr = &mesh->indexAttr();
    cellAttr = &mesh->attributes().attrMap<CCIndex, AuxinGradient::CellData>("AuxinCellData");
    edgeAttr = &mesh->attributes().attrMap<CCSignedIndex, AuxinGradient::EdgeData>("AuxinEdgeData");

    // Initialize the tissue
    if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw(QString("AuxinSolver::initialize Cannot make Tissue Process:" + parm("Tissue Process")));
    tissueProcess->initialize(parent);
    CCStructure &cs = tissue().cellStructure();

    // Initialize the solver derivatives provider
    if(!getProcess(parm("Growth Signal Derivatives"), auxinProcess))
      throw(QString("AuxinSolver::initialize Unable to make Growth Signal Derivatives:" + parm("Growth Signal Derivatives")));
    auxinProcess->initialize(tissue(), *indexAttr, *cellAttr, *edgeAttr);       
    for(CCIndex f : cs.faces()) {
      AuxinGradient::CellData &cCD = (*cellAttr)[f];
      if(cCD.fixedConc){
        mdxInfo << "the cell has fixed concentration equal to " << cCD.a << endl;
      }

    } 
    convergeThresh = parm("Converge Threshold").toDouble();
    // Register derivative providers and initialize
    clearDerivs();
    addDerivs(auxinProcess);
    initSolver(&tissue().dualGraph());

    return true;
  }

  bool AuxinSolver::step()
  {

    tissueProcess->tissue().updateGeometry();
    
    mesh->updateProperties(tissue().tissueName());
    mesh->updateProperties(tissue().tissueDualName());

    // RSS
    updateGeometry(tissueProcess->tissue().cellStructure(), mesh->indexAttr());
    mesh->updatePositions(tissueProcess->tissue().tissueDualName());   
    mesh->updatePositions(tissueProcess->tissue().tissueName());  
    // Update mesh points, edges, surfaces

    // Solve the system
    solve();

    // Copy the data to the attributes for visualization
    CCIndexDoubleAttr &auxinAttr = mesh->signalAttr<double>("Signal Auxin");
    CCSignedIndexDoubleAttr &pinAttr = currentMesh()->attributes().attrMap<CCSignedIndex, double>("Polarity");

    auxinProcess->maxA = 0;
    auxinProcess->minA = 0;
    double maxP = 0;
    CCStructure &cs = tissue().cellStructure();
    for(CCIndex f : cs.faces()) {
      AuxinGradient::CellData &cCD = (*cellAttr)[f];
      if(cCD.fixedConc){
        cCD.a = cCD.amount;
        mdxInfo << "the cell has fixed concentration equal to " << cCD.a  << endl;
      }
      if(auxinProcess->maxA < cCD.a)
        auxinProcess->maxA = cCD.a;
      if(auxinProcess->minA > cCD.a)
        auxinProcess->minA = cCD.a;
      auxinAttr[f] = cCD.a;
      pinAttr[f] = cCD.a;

      for(CCIndex e : cs.bounds(f)) {
        CCSignedIndex es(e, cs.ro(f, e));
        AuxinGradient::EdgeData &pED = (*edgeAttr)[es];

        if(maxP < pED.pin)
          maxP = pED.pin;

        pinAttr[es] = pED.pin;
      }
    }
    mdxInfo << "Max auxin: " << auxinProcess->maxA << " Max PIN:" << maxP << endl;

    mesh->updateProperties(tissue().tissueName());
    mesh->updateProperties(tissue().tissueDualName());

    // RSS
    updateGeometry(tissueProcess->tissue().cellStructure(), mesh->indexAttr());
    mesh->updatePositions(tissueProcess->tissue().tissueDualName());   
    mesh->updatePositions(tissueProcess->tissue().tissueName());        
    if(calcDXNorm() <= convergeThresh)
      return false;
    else
      return true;

  }

  // Return the Cell Tissue
  CellTissue &AuxinSolver::tissue()
  {
    if(!tissueProcess)
      throw(QString("AuxinGradient::tissue Tissue process empty"));
    return tissueProcess->tissue();
  }

  bool AuxinSetCellType::run()
  {
    cellType = AuxinGradient::stringToCellType(parm("Cell Type"));
    Amount = parm("Amount").toDouble();
    FixedConc = stringToBool(parm("Fixed Conc"));
    PropagateLaterally = stringToBool(parm("Propagate laterally"));

    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("AuxinSetCellType::run No current mesh"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("AuxinSetCellType::run Error, no cell complex selected"));

    CCStructure &cs = mesh->ccStructure(ccName);
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    AuxinGradient::CellDataAttr &cellAttr = mesh->attributes().attrMap<CCIndex,AuxinGradient::CellData>("AuxinCellData");

    for(CCIndex c : cs.faces()) {
      AuxinGradient::CellData &cCD = cellAttr[c];
      CCIndexData &cIdx = indexAttr[c];
      if(cIdx.selected) {
        cCD.cellType = cellType;
        cCD.amount = Amount;
        cCD.fixedConc = FixedConc;
        cCD.propagateLaterally = PropagateLaterally; 
      }
    }
    for(CCIndex c : cs.faces()) {
      AuxinGradient::CellData &cCD = cellAttr[c];
      CCIndexData &cIdx = indexAttr[c];
      mdxInfo << cCD.cellType << endl;

    }
    mesh->updateProperties(ccName);

    return true;
  }

  bool AuxinResetCellType::run()
  {
    cellType = AuxinGradient::GroundCell;
    Amount = 0;
    FixedConc = false;

    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("AuxinSetCellType::run No current mesh"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("AuxinSetCellType::run Error, no cell complex selected"));

    CCStructure &cs = mesh->ccStructure(ccName);
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    AuxinGradient::CellDataAttr &cellAttr = mesh->attributes().attrMap<CCIndex,AuxinGradient::CellData>("AuxinCellData");

    for(CCIndex c : cs.faces()) {
      AuxinGradient::CellData &cCD = cellAttr[c];
      CCIndexData &cIdx = indexAttr[c];
      //if(cIdx.selected) {
        cCD.cellType = cellType;
        cCD.amount = Amount;
        cCD.fixedConc = FixedConc;
      //}
    }

    for(CCIndex c : cs.faces()) {
      AuxinGradient::CellData &cCD = cellAttr[c];
      CCIndexData &cIdx = indexAttr[c];
      mdxInfo << cCD.cellType << endl;

    }
    mesh->updateProperties(ccName);

    return true;
  }


  void AuxinGradient::Subdivide::splitCellUpdate(Dimension dim, const CCStructure &cs, 
      const CCStructure::SplitStruct& ss, CCIndex otherP, CCIndex otherN,  double interpPos)
  {
    if(dim == 1) {
      if(!edgeAttr)
        std::cout << "AuxinGradient::Subdivide::splitCellUpdate Edge attribute not set" << std::endl;
      else
        (*edgeAttr)[ss.childP].pin = (*edgeAttr)[ss.childN].pin = (*edgeAttr)[ss.parent].pin;
    }

    if(dim == 2) {
      if(!cellAttr)
        std::cout << "AuxinGradient::Subdivide::splitCellUpdate Cell attribute not set" << std::endl;
      else {
    
         //CCIndexDataAttr &indexAttr = tissue->indexAttr();
         if(!indexAttr)
          std::cout << "AuxinGradient::Subdivide::splitCellUpdate Index data not set" << std::endl;

        AuxinGradient::CellData &pCD = (*cellAttr)[ss.childP];
        AuxinGradient::CellData &nCD = (*cellAttr)[ss.childN];
        AuxinGradient::CellData &rCD = (*cellAttr)[ss.parent];
        mdxInfo << "Prop Vector " << propagVector << endl;     
 
        pCD = nCD = rCD;
        if (rCD.fixedConc == true ){
          if(rCD.propFixConDivision == false)
          {
            mdxInfo << "I am in the propFixCondDivision false case " << endl;
            if ((rCD.posNearestFixed).norm() >= 1.e20)
            {
                 pCD.fixedConc = false; //the fixed father cell is alone. Keep one of the two daughters randomly
                 mdxInfo<< "isolated case of cells" << endl;
            }
            else
            {
                 /*CCIndexData &nIdx = (*indexAttr)[ss.childN];
                 CCIndexData &pIdx = (*indexAttr)[ss.childP];
                 Point3d p2n = nIdx.pos - pIdx.pos;
                 double aligned = p2n*(*propagVector);
                 if (aligned >= 1.e-6)
                   pCD.fixedConc = false;
                 else 
                   nCD.fixedConc = false; // the orthogonality is taken into account here */
                 CCIndexData &nIdx = (*indexAttr)[ss.childN];
                 CCIndexData &pIdx = (*indexAttr)[ss.childP];

                 if ((nIdx.pos- rCD.posNearestFixed).norm() < (pIdx.pos- rCD.posNearestFixed).norm()){
                   pCD.fixedConc = false;
                   mdxInfo<< "p cell set to false" << endl;

                 }
                 else  {
                   nCD.fixedConc = false;
                   mdxInfo<< "n cell set to false" << endl;
                 }
            }
            
            /*else
            {
               CCIndexData &nIdx = (*indexAttr)[ss.childN];
               CCIndexData &pIdx = (*indexAttr)[ss.childP];

               if ((nIdx.pos- rCD.posNearestFixed).norm() < (pIdx.pos- rCD.posNearestFixed).norm())
                 pCD.fixedConc = false;
               else  
                 nCD.fixedConc = false;
            }*/
           }
           else if (rCD.propagateLaterally == false){
             mdxInfo << "I am in the laterally false propagation case" << endl;
             //MassSpring::CellModelData  &rCMD = (*cellModelAttr)[ss.parent];
             CCIndexData &nIdx = (*indexAttr)[ss.childN];
             CCIndexData &pIdx = (*indexAttr)[ss.childP];
             Point3d joiningWall = (nIdx.pos - pIdx.pos)/(nIdx.pos - pIdx.pos).norm();
             mdxInfo << "joining Wall "<< joiningWall << endl;
             mdxInfo << "pol dir is " << rCD.copyPolDir << endl;
             if (fabs(joiningWall *rCD.copyPolDir) >= 0.5)
             {
               if ((nIdx.pos- rCD.posNearestFixed).norm() < (pIdx.pos- rCD.posNearestFixed).norm())
                 pCD.fixedConc = false;
               else  
                 nCD.fixedConc = false;
            }
          }

          }
        // What to do with PIN? For gradient it gets recalculated at every step.
      }
    }
  }	

  /*
   * Mass Spring
   */
  bool MassSpring::initialize(QWidget *parent)
  {
    // Get the current mesh
    mesh = currentMesh();
    if(!mesh)
      throw(QString("MassSpring::initialize No current mesh"));

    // Get the tissue, already initialized by the solver
    if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw(QString("MassSpring::initialize Cannot make tissue process") + parm("Tissue Process"));

    // Get the attribute maps
    indexAttr = &mesh->indexAttr();
    wallAttr = &mesh->attributes().attrMap<CCIndex,MassSpring::WallData>("MassSpringWallData");
    cellDAttr = &mesh->attributes().attrMap<CCIndex,MassSpring::CellModelData>("CellModelData");
    vtxDAttr  = &mesh->attributes().attrMap<CCIndex, MassSpring::VertexModelData>("VertexModelData");

    delta = parm("Delta").toDouble();
    springK = parm("Spring Constant").toDouble();
    pressure = parm("Turgor Pressure").toDouble();
    mdxInfo << "The pressure is: " << pressure << endl;
    // assign global or individual cell based pressure and stiffness values

    return true;
  }
  
  bool SetDirichlet::initialize(QWidget *parent)
  {
    dirichletFixed = stringToPoint3u(parm("Dirichlet"));
    nodalNeumann = stringToPoint3d(parm("Nodal Neumann"));

    return true;
  }

  bool SetDirichlet::run()
  {
    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("MassSpring::run No current mesh"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("SetDirichlet::run Error, no cell complex selected"));

    CCStructure &cs = mesh->ccStructure(ccName);
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    MassSpring::VertexModelAttr  &vMAttr =  mesh->attributes().attrMap<CCIndex,MassSpring::VertexModelData>("VertexModelData");
    dirichletAttr = &mesh->attributes().attrMap<CCIndex, fem::Dirichlet>(parm("Dirichlet Attribute"));
	

    forall(CCIndex c, cs.cellsOfDimension(0)) {
      CCIndexData &cIdx = indexAttr[c];
      if(cIdx.selected) {
        (vMAttr)[c].dirichlet = dirichletFixed;
        (vMAttr)[c].Neumann = nodalNeumann;
        (*dirichletAttr)[c].labels = dirichletFixed;
      }
    }
    mesh->updateProperties(ccName);

    return true;
  }

  
  bool SpringSetCellType::run()
  {
    cellTypeSpr = MassSpring::stringToCellTypeSpr(parm("Cell Type"));
    maxArea = parm("maxArea").toDouble();
    isoGrowth = stringToBool(parm("Isotropic Growth"));
    paralGrowth = parm("Parallel Growth to Field").toDouble();
    orthoGrowth = parm("Orthogonal Growth to Field").toDouble();
    cellPressure = parm("Cell Pressure").toDouble();
    //cellStiffness = parm("Cell Stiffness").toDouble();
    strainThresh = parm("Strain threshold").toDouble();
    strBasedGrRate = parm("Strain based growth rate").toDouble();
    specialGrCon = parm("Special growth rate intensity").toDouble();
    growthRateStrParall = parm("Strain based growth rate par to field").toDouble();
    growthRateStrOrtho = parm("Strain based growth rate ortho to field").toDouble();


    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("SpringSetCellType::run No current mesh"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("SpringSetCellType::run Error, no cell complex selected"));

    CCStructure &cs = mesh->ccStructure(ccName);
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    MassSpring::CellModelAttr &cellAttr = mesh->attributes().attrMap<CCIndex,MassSpring::CellModelData>("CellModelData");

    for(CCIndex c : cs.faces()) {
      MassSpring::CellModelData &cCD = cellAttr[c];
      CCIndexData &cIdx = indexAttr[c];
      if(cIdx.selected) {
        cCD.cellTypeSpr = cellTypeSpr;
        cCD.maxAreaDiv = maxArea;
        cCD.cellIsoGrowth =isoGrowth;
        cCD.cellParallGrowth = paralGrowth;
        cCD.cellOrthoGrowth = orthoGrowth;
        cCD.cellPressure = cellPressure;
        //cCD.cellStiffness = cellStiffness;
        cCD.strainThreshold = strainThresh;
        cCD.strBasedGrRate = strBasedGrRate;
        cCD.specialGrCon = specialGrCon;
        cCD.growthRateStrParall = growthRateStrParall;
        cCD.growthRateStrOrtho  = growthRateStrOrtho;
      }
    }
    mesh->updateProperties(ccName);

    return true;
  }
  bool OnlyPericlinalGrowthL1::run()
  {
      onlyPericlinalGrL1 = stringToBool(parm("Only Periclinal Growth in L1"));
      if (!onlyPericlinalGrL1)
          return true;
      else{
        Mesh *mesh = currentMesh();
        if(!mesh)
          throw(QString("OnlyPericlinalGrowthL1::run No current mesh"));
     
        QString ccName = mesh->ccName();
        if(ccName.isEmpty())
          throw(QString("OnlyPericlinalGrowthL1::run Error, no cell complex selected"));
     

        MassSpring::WallDataAttr  &wDAttr =  mesh->attributes().attrMap<CCIndex,MassSpring::WallData>("MassSpringWallData");
        MassSpring::CellModelAttr &cellDAttr = mesh->attributes().attrMap<CCIndex,MassSpring::CellModelData>("CellModelData");

  
        CCStructure &cs = mesh->ccStructure(ccName);
        CCIndexDataAttr &indexAttr = mesh->indexAttr();
        for(CCIndex c : cs.faces()) {
          MassSpring::CellModelData &cCD = cellDAttr[c];
          CCIndexData &cIdx = indexAttr[c];
          if (cCD.cellTypeSpr == MassSpring::L1)
              for (CCIndex edge : cs.incidentCells(c,1)){
                if (cs.cobounds(edge).size() == 2){
                    uint boundaryL1Cells = 0;
                    for (CCIndex coBoundCell : cs.cobounds(edge))
                       if ((cellDAttr)[coBoundCell].cellTypeSpr == MassSpring::L1)
                         boundaryL1Cells++;
                    if(boundaryL1Cells == 2)
                        (wDAttr)[edge].noGrowth = true;
                }        
              }
            
        }         
        return true; 
     }   
  }
    
  bool SpringResetCellType::run()
  {
    cellTypeSpr = MassSpring::NormalCell;

    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("SpringReetCellType::run No current mesh"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("SpringReetCellType::run Error, no cell complex selected"));

    CCStructure &cs = mesh->ccStructure(ccName);
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    MassSpring::CellModelAttr &cellAttr = mesh->attributes().attrMap<CCIndex,MassSpring::CellModelData>("CellModelData");

    for(CCIndex c : cs.faces()) {
      MassSpring::CellModelData &cCD = cellAttr[c];
      CCIndexData &cIdx = indexAttr[c];
      //if(cIdx.selected) {
        cCD.cellTypeSpr = cellTypeSpr;
      //}
    }
    mesh->updateProperties(ccName);

    return true;
  }


  bool SetWallStiffness::initialize(QWidget *parent)
  {
    specialStiffness = parm("Wall Stiffness").toDouble();
    stiffnessBasedMorpho = stringToBool(parm("Wall Stiffness based on morphogen"));
    maxStiffness = parm("Wall Stiffness Max").toDouble();
    minStiffness = parm("Wall Stiffness Min").toDouble();
    proportionality = parm("Proportionality type");

    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("Set Wall Stiffness::run No current mesh"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("SetWallStiffness::run Error, no cell complex selected"));



    /*if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw(QString("SetWallStiffness::initialize Cannot make Tissue Process:" + parm("Tissue Process")));
    tissueProcess->initialize(parent);*/

    //CCStructure &cs = mesh->ccStructure(ccName);
    //CCIndexDataAttr &indexAttr = mesh->indexAttr();
    //MassSpring::VertexModelAttr  &vMAttr =  mesh->attributes().attrMap<CCIndex,MassSpring::VertexModelData>("VertexModelData");
    //MassSpring::WallDataAttr  &wDAttr =  mesh->attributes().attrMap<CCIndex,MassSpring::WallData>("MassSpringWallData");
    //AuxinGradient::CellDataAttr &cellAttr = mesh->attributes().attrMap<CCIndex,AuxinGradient::CellData>("AuxinCellData");

// Initialize the solver derivatives provider
    if(!getProcess(parm("Growth Signal Gradient"), auxinGradientStiff))
      throw(QString("SetWallStiffness::run Unable to retrieve the Auxin Gradient process"));

   return true;
  }
 
  bool SetWallStiffness::run()
  {
    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("Set Wall Stiffness::run No current mesh"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("SetWallStiffness::run Error, no cell complex selected"));

    CCStructure &cs = mesh->ccStructure(ccName);
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    MassSpring::VertexModelAttr  &vMAttr =  mesh->attributes().attrMap<CCIndex,MassSpring::VertexModelData>("VertexModelData");
    MassSpring::WallDataAttr  &wDAttr =  mesh->attributes().attrMap<CCIndex,MassSpring::WallData>("MassSpringWallData");
    AuxinGradient::CellDataAttr &cellAttr = mesh->attributes().attrMap<CCIndex,AuxinGradient::CellData>("AuxinCellData");


    for(CCIndex e : cs.edges()) {
         CCIndexData &eIdx = indexAttr[e];
         MassSpring::WallData &eD = wDAttr[e];
         for(CCIndex f : cs.cobounds(e)){ 
           if (indexAttr[f].selected == true) //if there is a selection a special stiffness will be assigned
           { 
              if ((wDAttr)[e].wallStiffnessAssigned == true)
                (wDAttr)[e].wallStiffness = ((wDAttr)[e].wallStiffness + specialStiffness)/2.; //an edge in 2D can have at most 2 cobound cells 
              else {
                (wDAttr)[e].wallStiffnessAssigned = true;
                (wDAttr)[e].wallStiffness = specialStiffness; // not doing averaging with non-selected walls

              }
             }

             if (stiffnessBasedMorpho == true)
             {
                 AuxinGradient::CellData &aCD = cellAttr[f];
                 double stiffTemp = 0;
                 if (proportionality == "Positive")
                   stiffTemp = minStiffness + ((aCD.a - auxinGradientStiff->minA) * (maxStiffness-minStiffness)/(auxinGradientStiff->maxA-auxinGradientStiff->minA));
                 else
                   stiffTemp = minStiffness + ((-aCD.a + auxinGradientStiff->maxA) * (maxStiffness-minStiffness)/(auxinGradientStiff->maxA-auxinGradientStiff->minA));
                 if ((wDAttr)[e].wallStiffnessAssigned == true)
                   (wDAttr)[e].wallStiffness = ((wDAttr)[e].wallStiffness + stiffTemp)/2.; //an edge in 2D can have at most 2 cobound cells     
                 else {
                   (wDAttr)[e].wallStiffnessAssigned = true;
                   (wDAttr)[e].wallStiffness = stiffTemp; // not doing averaging with non-selected walls
                 }
             }
                   
           }
            mdxInfo << "The stiffness assigned is: " << (wDAttr)[e].wallStiffness << endl;
         } 
              
    
    mesh->updateProperties(ccName);
    return true;
  }

  bool ResetWallStiffness::run()
  {

    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("Set Wall Stiffness::run No current mesh"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("SetWallStiffness::run Error, no cell complex selected"));

    CCStructure &cs = mesh->ccStructure(ccName);
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    MassSpring::WallDataAttr  &wDAttr =  mesh->attributes().attrMap<CCIndex,MassSpring::WallData>("MassSpringWallData");
   
    forall(CCIndex e, cs.cellsOfDimension(1)){
      CCIndexData &eIdx = indexAttr[e];
      
      {
          (wDAttr)[e].wallStiffnessAssigned = false;
          (wDAttr)[e].wallStiffness = 0;
          mdxInfo << "I have reset all the stiffnesses to 0, no wallStiffnessAssigned equal to true" << endl;
      }
    }
   
    mesh->updateProperties(ccName);
    return true;
  }

 
  bool PolarizerSetFixedConc::run()
  {
    fixedConcentration = parm("Concentration").toDouble();


    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("PolarizerSetFixedConc::run No current mesh"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("PolarizerSetFixedConc::run Error, no cell complex selected"));

    CCStructure &cs = mesh->ccStructure(ccName);
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    MassSpring::CellModelAttr &cellDAttr = mesh->attributes().attrMap<CCIndex,MassSpring::CellModelData>("CellModelData");

    for(CCIndex c : cs.faces()) {
      MassSpring::CellModelData &cCD = cellDAttr[c];
      CCIndexData &cIdx = indexAttr[c];
      if(cIdx.selected) {
        cCD.fixDiffAnisoConc = true;
        cCD.fixedConcVal = fixedConcentration;
      }
    }
    for(CCIndex c : cs.faces()) {
      MassSpring::CellModelData &cCD = cellDAttr[c];
      CCIndexData &cIdx = indexAttr[c];
      mdxInfo << cCD.fixedConcVal << endl;
    }

    mesh->updateProperties(ccName);

    return true;
  }

  bool PolarizerResetFixedConc::run()
  {

    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("PolarizerResetFixedConc::run No current mesh"));

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw(QString("PolarizerResetFixedConc::run Error, no cell complex selected"));

    CCStructure &cs = mesh->ccStructure(ccName);
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    MassSpring::CellModelAttr &cellDAttr = mesh->attributes().attrMap<CCIndex,MassSpring::CellModelData>("CellModelData");

    for(CCIndex c : cs.faces()) {
      MassSpring::CellModelData &cCD = cellDAttr[c];
      CCIndexData &cIdx = indexAttr[c];
      //if(cIdx.selected) {
        cCD.fixDiffAnisoConc = false;
        cCD.fixedConcVal = 0.;
      //}
    }
    for(CCIndex c : cs.faces()) {
      MassSpring::CellModelData &cCD = cellDAttr[c];
      CCIndexData &cIdx = indexAttr[c];
      mdxInfo << cCD.fixedConcVal << endl;
    }

    mesh->updateProperties(ccName);

    return true;
  }
 

  bool PolarizerDiffusion::initialize(QWidget *parent) 
  {
    mesh = currentMesh();
    if(!mesh)
      throw(QString("PolarizerDiffusion::initialize No current mesh"));

    // Get the tissue, already initialized by the solver
    if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw(QString("Polarizer Diffusion::initialize Cannot make tissue process") + parm("Tissue Process"));

    // Get the attribute maps
    indexAttr = &mesh->indexAttr();
    //wallAttr = &mesh->attributes().attrMap<CCIndex,MassSpring::WallData>("MassSpringWallData");
    cellDAttr = &mesh->attributes().attrMap<CCIndex,MassSpring::CellModelData>("CellModelData");

    DiffConst = parm("Diffusion Constant").toDouble();
    return true;
  }
  // polarizertropy initialize derivatives
  /*void PolarizerDiffusion::initDerivs(SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr)
  {
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    //WallDataAttr &wallAttr = mesh->attributes().attrMap<CCIndex,WallData>("MassSpringWallData");
    //CCStructure &cs = tissueProcess->tissue().cellStructure();
    MassSpring::VertexModelAttr &vMAttr = mesh->attributes().attrMap<CCIndex,MassSpring::VertexModelData>("VertexModelData");
  }*/
  // Calculate derivative for the solver
  void PolarizerDiffusion::calcDerivatives(const SolverT &solver, VertexAttr &vertexAttr)
  {
    CCStructure &dg = tissueProcess->tissue().dualGraph();
    //CCIndexDataAttr &indexAttr = mesh->indexAttr();
    //MassSpring::VertexModelAttr &vMAttr = mesh->attributes().attrMap<CCIndex,MassSpring::VertexModelData>("VertexModelData");
    //CCStructure &cs = tissueProcess->tissue->cellStructure();
    //CCStructure &dg = tissueProcess->tissue->dualGraph();

    // First clear derivatives
    for(CCIndex c : dg.vertices())
      vertexAttr[c].dx.zero();

    // Now calculate fixed concentration
    for(CCIndex c : dg.vertices()) {
      MassSpring::CellModelData &cCD = (*cellDAttr)[c];
      if(cCD.fixDiffAnisoConc)
        cCD.concentration = cCD.fixedConcVal;
    }
    // Now calculate
    for(CCIndex c : dg.vertices()) {
      MassSpring::CellModelData &cCD = (*cellDAttr)[c];
      CCIndexData &cIdx = (*indexAttr)[c];
      SolverT::SolverVertexData &cSD = vertexAttr[c];

      // Calculate diffusion
      for(const Flip &flip : neighborVertexFlips(dg, c)) {
        CCIndex w = flip.interior;
        CCIndex n = flip.otherFacet(c);
  
        // Diffusion
        MassSpring::CellModelData &nCD = (*cellDAttr)[n];
        //SolverT::SolverVertexData &nSD = vertexAttr[n];
        CCIndexData &wIdx = (*indexAttr)[w];
        if (cCD.fixDiffAnisoConc == true)
          cSD.dx[0] = 0.;
        else{ 
         cSD.dx[0] += DiffConst * (nCD.concentration - cCD.concentration) * wIdx.measure / cIdx.measure;

        //mdxInfo<< "cIdx measure " <<  cIdx.measure  << endl;
        //mdxInfo<< "cIdx position " << cIdx.pos << endl;
        //mdxInfo<< "cSD dx " <<   cSD.dx[0] << endl;
        }
      }
    }
  }

  // Return the Cell Tissue
  /*CellTissue &PolarizerDiffusionSolver::tissue()
  {
    if(!tissueProcess)
      throw(QString("PolarizerDiffusionSolver::tissue Tissue process empty"));
    return tissueProcess->tissue();
  }*/

  bool PolarizerDiffusionSolver::initialize(QWidget *parent)
  {
    mesh = currentMesh();
    if(!mesh)
      throw(QString("AuxinSolver::initialize No current mesh"));

    indexAttr = &mesh->indexAttr();
    cellDAttr = &mesh->attributes().attrMap<CCIndex, MassSpring::CellModelData>("CellModelData");
    //wallAttr = &mesh->attributes().attrMap<CCSignedIndex, MassSpring::WallData>("MassSpringWallData");

    // Initialize the tissue
    if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw(QString("PolarizerDiffusion::initialize Cannot make Tissue Process:" + parm("Tissue Process")));
    tissueProcess->initialize(parent);

    // Initialize the solver derivatives provider
    if(!getProcess(parm("Polarizer Diffusion Process"), polarizerProcess))
      throw(QString("PolarizerSolver::initialize Unable to make Polarizer Process:" + parm("Polarizer Process")));
    //polarizertropyProcess->initialize(tissue(), *indexAttr, *cellDAttr, *wallAttr);
    polarizerProcess->initialize(parent);        



    convergeThresh = parm("Converge Threshold").toDouble();
    // Register derivative providers and initialize
    clearDerivs();
    addDerivs(polarizerProcess);
    initSolver(&tissueProcess->tissue().dualGraph());	


    return true;
  }

  bool PolarizerDiffusionSolver::step()
  { 


    tissueProcess->tissue().updateGeometry();
    
    mesh->updateProperties(tissueProcess->tissue().tissueName());
    mesh->updateProperties(tissueProcess->tissue().tissueDualName());

    // RSS
    updateGeometry(tissueProcess->tissue().cellStructure(), mesh->indexAttr());
    mesh->updatePositions(tissueProcess->tissue().tissueDualName());   
    mesh->updatePositions(tissueProcess->tissue().tissueName()); 
    // Solve the system
    solve();

    // Copy the data to the attributes for visualization
    //CCIndexDoubleAttr &polarizerAttr = mesh->signalAttr<double>("Signal");
    //CCSignedIndexDoubleAttr &pinAttr = currentMesh()->attributes().attrMap<CCSignedIndex, double>("Polarity");

    double maxC = 0;
    //double maxP = 0;
    CCStructure &cs = tissueProcess->tissue().cellStructure();
    for(CCIndex f : cs.faces()) {
      MassSpring::CellModelData &cCD = (*cellDAttr)[f];
      if(cCD.fixDiffAnisoConc){
        cCD.concentration = cCD.fixedConcVal;
        //mdxInfo << "the cell has fixed concentration" << endl;
      }
      if(maxC < cCD.concentration)
        maxC = cCD.concentration;
    //  polarizerAttr[f] = cCD.concentration;

    }
    mdxInfo << "Max concentration polarizer: " << maxC << endl;
    

    mesh->updateProperties(tissueProcess->tissue().tissueName());
    mesh->updateProperties(tissueProcess->tissue().tissueDualName());

    // RSS
    updateGeometry(tissueProcess->tissue().cellStructure(), mesh->indexAttr());
    mesh->updatePositions(tissueProcess->tissue().tissueDualName());   
    mesh->updatePositions(tissueProcess->tissue().tissueName()); 
    if(calcDXNorm() <= convergeThresh)
    {
      //compute the pseudo-gradient of this diffusion, the vector used to set polarizer
      for(CCIndex f : cs.faces()) {
       Point3d gradientConc = Point3d(0., 0., 0.);
       MassSpring::CellModelData &cCD = (*cellDAttr)[f];
       for(CCIndex n : cs.neighbors(f)) {
         MassSpring::CellModelData &nCD = (*cellDAttr)[n];
         Point3d distanceVector = ((*indexAttr)[f].pos - (*indexAttr)[n].pos)/((*indexAttr)[f].pos - (*indexAttr)[n].pos).norm();
         gradientConc += distanceVector*(cCD.concentration - nCD.concentration);///((*indexAttr)[f].pos.x() - (*indexAttr)[n].pos.x());
       }
       if (gradientConc.norm() >= 1.e-12){
         gradientConc *= 1./gradientConc.norm();
         cCD.polarizerDir = gradientConc;
       }
       else
         cCD.polarizerDir = 0.5*Point3d(sqrtf(2.), sqrtf(2.), 0.);
      }
      return false;

    }
    else
      return true;

  }


  // use it to compute the pressure attribute
  Point3d MassSpring::calcDerivsCell(const CCStructure &cs, const CCIndexDataAttr &indexAttr, CCIndex f, CCIndex e, CCIndex v)
  {
    Point3d dx;
    Point3d vPos = indexAttr[v].pos;

    //double pressure = 0.3;
    // Find the centriod
    Point3d centroid = indexAttr[f].pos;
    Point3d edgeMid;
    auto eb = cs.edgeBounds(e);
    Point3d nPos;
    if(eb.first == v)
      nPos = indexAttr[eb.second].pos;
    else
      nPos = indexAttr[eb.first].pos;
    edgeMid = (nPos + (vPos - nPos)/2.);
      
    Point3d pressureVersor = (edgeMid - centroid);
    Point3d pressureDir = Point3d(0., 0., 1.)^(vPos - nPos);
    pressureDir *= 1./norm(pressureDir);
    if (pressureDir*pressureVersor < 0.)
        pressureDir *= -1.;
     
    MassSpring::CellModelData &cCD = (*cellDAttr)[f];
    if (cCD.cellTypeSpr == NormalCell){ 
     dx += pressure * norm(vPos - nPos) * pressureDir; 
     cCD.cellPressure = pressure;
    }
    else
     dx += cCD.cellPressure * norm(vPos - nPos) * pressureDir; 
    //mdxInfo << "the pressure is: " << pressure << endl;
    return dx;
  }

  // this function computes the force acting on the vtx v as member of the edge e. Its connected vtx n 
  Point3d MassSpring::calcDerivsWall(const CCStructure &cs, const CCIndexDataAttr &indexAttr, CCIndex e, CCIndex v)
  {
    Point3d dx;
    //double springK = 100.; // to become a model parameter
    Point3d vPos = indexAttr[v].pos;
    auto eb = cs.edgeBounds(e);
    Point3d nPos;
    if(eb.first == v)
      nPos = indexAttr[eb.second].pos;
    else
      nPos = indexAttr[eb.first].pos;
    const WallData &eD = (*wallAttr)[e];

    double currLength =  norm(vPos - nPos);
    double localSpringK = springK;
 
    if (eD.wallStiffnessAssigned == true)
        localSpringK = eD.wallStiffness;
    //mdxInfo << "The local spring is : " << localSpringK << endl;
    dx += localSpringK * ((eD.restLength - currLength)/eD.restLength) *(vPos-nPos)/(vPos-nPos).norm();
    //mdxInfo << "the spring K is " << springK << endl;
    return dx;
  }

  void MassSpring::calcDerivatives(const SolverT &solver, CCIndex v, VectorT &dx)
  {
    const CCStructure &cs = tissueProcess->tissue().cellStructure();
    CCIndexUSet edgeSet;
    for(auto flip : cs.matchV(v, CCIndex::Q, CCIndex::Q, CCIndex::Q)) {
      CCIndex f = flip.interior;
      if(f.isPseudocell())
        continue;

      // Save attached edges
      edgeSet.insert(flip.facet[0]);
      edgeSet.insert(flip.facet[1]);
      // Calculate the derivatives contributions from
      dx += calcDerivsCell(cs, *indexAttr, f, flip.facet[0], v);
      dx += calcDerivsCell(cs, *indexAttr, f, flip.facet[1], v);
    }
    for(CCIndex e : edgeSet)
      dx += calcDerivsWall(cs, *indexAttr, e, v);
    MassSpring::VertexModelData &vMD = (*vtxDAttr)[v];
    //I think Neumann is overridden by Dirichlet
    if ((fabs(vMD.Neumann.x()) >= (1.e-3)) or (fabs(vMD.Neumann.y())>= (1.e-3))  or (fabs(vMD.Neumann.z())>= (1.e-3)) )
       dx += vMD.Neumann; 


    return;
  }

  void MassSpring::initDerivs(SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr)
  {
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    WallDataAttr &wallAttr = mesh->attributes().attrMap<CCIndex,WallData>("MassSpringWallData");
    CCStructure &cs = tissueProcess->tissue().cellStructure();

    // Fill in empty target areas and face vertex position pointers
//    for(CCIndex f : cs.faces()) {
//      
//      double area = indexAttr[f].measure;
//      if(fD.targetArea <= 0)
//        fD.targetArea = area;
//
//      // Setup the pointers to positions for area calculation later
//      fD.posVec.clear();
//      fD.vtxVec.clear();
//      CellTuple faceVtx(cs,f);
//      CCIndex firstVertex = faceVtx[0];
//      do {
//        fD.posVec.push_back(&indexAttr[faceVtx[0]].pos);
//        fD.vtxVec.push_back(faceVtx[0]);
//        faceVtx.flip(0,1);
//      } while(faceVtx[0] != firstVertex);
//    } 
    //MassSpring::VertexModelAttr &vMAttr = mesh->attributes().attrMap<CCIndex,MassSpring::VertexModelData>("VertexModelData");
    //set the Dirichlet boundary 
    for(CCIndex v : cs.cellsOfDimension(0)) {
      MassSpring::VertexModelData &vMD = (*vtxDAttr)[v];
      if (vMD.dirichlet.x() != 0)
       vAttr[v].dirichlet.x() = 1;
      if (vMD.dirichlet.y() != 0)
       vAttr[v].dirichlet.y() = 1;
      if (vMD.dirichlet.z() != 0)
       vAttr[v].dirichlet.z() = 1;
    }
     
    for(CCIndex e : cs.edges()) {
      WallData &eD = wallAttr[e];

      // Set rest length if empty
      if(eD.restLength <= 0) {
        auto eb = cs.edgeBounds(e);
        double length = norm(indexAttr[eb.first].pos - indexAttr[eb.second].pos);
        //if(wallPreStress)
	      //eD.restLength = length * (1.0 - wallThresh);
        //else
	eD.restLength = length;
      }
     //mdxInfo<< "Rest Ln " << eD.restLength << endl;
    }
  }

  bool MassSpringSolver::initialize(QWidget *parent)
  {
    // Get the current mesh
    mesh = currentMesh();
    if(!mesh)
      throw(QString("MassSpringSolver::initialize No current mesh"));

    // Initialize the tissue
    if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw(QString("MassSpringSolver::initialize Cannot make Tissue Process:") + parm("Tissue Process"));
    tissueProcess->initialize(parent);

    // Initialize the solver derivs process
    if(!getProcess(parm("Mass Spring Process"), massSpringProcess))
      throw(QString("MassSpringGrowth::initialize Cannot make Mass Spring Process:") + parm("Mass Spring Process"));
    massSpringProcess->initialize(parent);

    clearDerivs();
    addDerivs(massSpringProcess);
    initSolver(&tissueProcess->tissue().cellStructure());	

    // Get the parms
    convergeThresh = parm("Converge Threshold").toDouble();
    mdxInfo << "The mass spring Solver has been initialized" << endl;
    return true;
  }

  bool MassSpringSolver::step()
  {
    if(!mesh)
      throw(QString("MassSpringSolver::step No current mesh"));

    if(!tissueProcess)
      throw(QString("MassSpringSolver::step No Tissue process"));

    if(!massSpringProcess)
      throw(QString("MassSpringSolver::step No MassSPring process"));

    solve();
    // RSS
    updateGeometry(tissueProcess->tissue().cellStructure(), mesh->indexAttr());
    mesh->updatePositions(tissueProcess->tissue().tissueName());         
    //mdxInfo << "I am in mass spring solver" << endl;
    if(calcDXNorm() <= convergeThresh){
       mdxInfo << " I detect mass spring convergence" << endl;
//      mdxInfo << "I see I have converged " << endl;
//      //save a snapshot of the simulation
//      static int screenShotCount = 0;
//      // Take screenshot
//      QString fileName = QString("OvuleGrowth-%1.JPG").arg(screenShotCount++, 4, 10, QChar('0'));
//      takeSnapshot(fileName);
//      mdxInfo << "Snapshot saved to file:" << fileName << endl;
      return false;
    }
    else
      return true;
  }

  bool MassSpringGrowth::initialize(QWidget *parent)
  {
    // Get the current mesh
    mesh = currentMesh();
    if(!mesh)
      throw(QString("MassSpringGrowth::initialize No current mesh"));

    // Grab the tissue process
    if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw(QString("MassSpringGrowth::initialize Cannot make tissue process:") + parm("Tissue Process"));
    tissueProcess->initialize(parent); // Is this required? How to avoid multiple initializations

    growthDt = parm("Growth Signal-Based Dt").toDouble();
    growthRate = parm("Growth Rate Strain").toDouble();
    strainThresh = parm("Strain Threshold").toDouble();
    isoGrowth = stringToBool(parm("Isotropic Growth"));
    paralGrowth = parm("Parallel Growth to Field").toDouble();
    orthoGrowth = parm("Orthogonal Growth to Field").toDouble();
    growthRateStrParall = parm("Strain based growth rate par to field").toDouble();
    growthRateStrOrtho = parm("Strain based growth rate ortho to field").toDouble();

    return true;
  }

  bool MassSpringGrowth::step()
  {
    if(!mesh)
      throw(QString("MassSpringSolver::step No current mesh"));

    if(!tissueProcess)
      throw(QString("MassSpringSolver::step No Tissue process"));

    strainThresh = parm("Strain Threshold").toDouble();
    CCStructure &cs = tissueProcess->tissue().cellStructure();

    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    //MassSpring::CellDataAttr &cellSpringDAttr = mesh->attributes().attrMap<CCIndex,MassSpring::>("MassSpringWallData");
    MassSpring::CellModelAttr &cellSpringAttr = mesh->attributes().attrMap<CCIndex, MassSpring::CellModelData>("CellModelData");
    MassSpring::WallDataAttr &wallAttr = mesh->attributes().attrMap<CCIndex,MassSpring::WallData>("MassSpringWallData");
    AuxinGradient::CellDataAttr &cellADAttr = mesh->attributes().attrMap<CCIndex,AuxinGradient::CellData>("AuxinCellData");
    
    for(CCIndex e : cs.edges()) {
        MassSpring::WallData &eD = wallAttr[e];
        if (eD.noGrowth == true)
          continue;
        eD.wallGrowthFactor = 0;
        //double growthRateAver = 0;
        //double strainThreshAver = 0.;
        eD.wallStrainGrowthFactor = 0.;
        eD.wallStrainThreshold = 0.;
        eD.wallStrainParallelGrowthFactor = 0.;
        eD.wallStrainOrthoGrowthFactor = 0.;
        int numIncidFaces = 0;
        for(CCIndex face: cs.incidentCells(e,2))
        {
          ++numIncidFaces;
          //CellData &cD = auxinGridAttr[face];
          MassSpring::CellModelData &cCD = cellSpringAttr [face];
          
          if (cCD.cellTypeSpr == MassSpring::NormalCell)
          {
             eD.wallStrainGrowthFactor += growthRate;
             eD.wallStrainThreshold += strainThresh;
             eD.wallGrowthFactor += (cellADAttr)[face].a;
             cCD.strBasedGrRate = growthRate;
             cCD.strainThreshold = strainThresh;
          } 
          else
          {
             //mdxInfo << "I am dealing with a special cell" << endl;
             eD.wallStrainGrowthFactor += cCD.strBasedGrRate;
             eD.wallStrainThreshold += cCD.strainThreshold;
             if (cCD.specialGrCon > 0.)
              eD.wallGrowthFactor += cCD.specialGrCon;
             else 
              eD.wallGrowthFactor += (cellADAttr)[face].a;
          }
        }
        eD.wallGrowthFactor *=1./numIncidFaces;
        eD.wallStrainGrowthFactor *= 1./numIncidFaces;
        eD.wallStrainThreshold *= 1./numIncidFaces;
   
        auto eb = cs.edgeBounds(e);
        double length = norm(indexAttr[eb.first].pos - indexAttr[eb.second].pos);
  	double strain = (length - eD.restLength)/eD.restLength;   
  	if(strain > strainThresh){
           eD.restLength += eD.restLength*(strain-eD.wallStrainThreshold) *  eD.wallStrainGrowthFactor;
           if(eD.restLength > length) {
            mdxInfo << "Restlength greater than length, rate too high?" << endl;
            eD.restLength = length;
           }

        }
        if (isoGrowth)
         eD.restLength +=  eD.restLength*growthDt*eD.wallGrowthFactor;
        else
        {
         //do polarizertropic growth - polarization field -based    
         //grow in each cell  
         Point3d polarizerDir = Point3d(0., 0., 0.);
         double averParaGrowth = 0.;
         double averOrthoGrowth = 0.;
         double averParaGrowthStrain = 0.;
         double averOrthoGrowthStrain = 0.;

         for(CCIndex face: cs.incidentCells(e,2))
         {
          MassSpring::CellModelData &cSD = cellSpringAttr[face];
          polarizerDir += cSD.polarizerDir; // the polarizer direction is already normalized
          if (cSD.cellTypeSpr != MassSpring::NormalCell)

          {
            averParaGrowth += cSD.cellParallGrowth;
            averOrthoGrowth += cSD.cellOrthoGrowth;
            averParaGrowthStrain += cSD.growthRateStrParall;
            averOrthoGrowthStrain += cSD.growthRateStrOrtho;

          }
          else
          {
            averParaGrowth += paralGrowth;
            averOrthoGrowth += orthoGrowth;
            averParaGrowthStrain += growthRateStrParall;
            averOrthoGrowthStrain += growthRateStrOrtho;
            cSD.cellParallGrowth = paralGrowth;
            cSD.cellOrthoGrowth = orthoGrowth;
            cSD.growthRateStrParall = growthRateStrParall;
            cSD.growthRateStrOrtho = growthRateStrOrtho;

          }
         }
         averParaGrowth *= 1./numIncidFaces;
         averOrthoGrowth *= 1./numIncidFaces;
         averParaGrowthStrain *= 1./numIncidFaces;
         averOrthoGrowthStrain *= 1./numIncidFaces;

         double polarizerDirNorm = norm(polarizerDir);
         polarizerDir *= 1./polarizerDirNorm;
         Point3d currentVersor = (indexAttr[eb.first].pos - indexAttr[eb.second].pos)/norm(indexAttr[eb.first].pos - indexAttr[eb.second].pos);
         // get the projection of the versor along the parallel direction
         double paralProj =  fabs(currentVersor*polarizerDir);
         // get the projection of the versor along the orthogonal direction
         // I do not care about the orientation of this vector, only its direction
         Point3d polarizerOrthDir = polarizerDir ^ Point3d(0., 0., 1.);
         polarizerOrthDir *= 1./norm(polarizerDir ^ Point3d(0., 0., 1.));
         double orthoProj = fabs(currentVersor*polarizerOrthDir);
         if (paralProj/orthoProj >=0.6) 
             eD.restLength +=  (eD.restLength*averParaGrowth*growthDt*eD.wallGrowthFactor) + (eD.restLength*(strain-eD.wallStrainThreshold) * averParaGrowthStrain);
         else if (paralProj/orthoProj <=0.4)
             eD.restLength +=  (eD.restLength*averOrthoGrowth*growthDt*eD.wallGrowthFactor)+ (eD.restLength*(strain-eD.wallStrainThreshold) * averOrthoGrowthStrain);

         else 
          eD.restLength +=  0.5 * ((eD.restLength*paralProj*averParaGrowth*growthDt*eD.wallGrowthFactor) + (eD.restLength*orthoProj*averOrthoGrowth*growthDt*eD.wallGrowthFactor) + (eD.restLength*(strain-eD.wallStrainThreshold) * averParaGrowthStrain) + (eD.restLength*(strain-eD.wallStrainThreshold) * averOrthoGrowthStrain));
 
       }
       //mdxInfo << "restLength after growth is: " << eD.restLength << endl; 
      };

    return true;
  }

 
 /* bool MassSpringGrowth::step()
  {
    if(!mesh)
      throw(QString("MassSpringSolver::step No current mesh"));

    if(!tissueProcess)
      throw(QString("MassSpringSolver::step No Tissue process"));

    //strainThresh = parm("Strain Threshold").toDouble();
    CCStructure &cs = tissueProcess->tissue().cellStructure();

    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    MassSpring::WallDataAttr &wallAttr = mesh->attributes().attrMap<CCIndex,MassSpring::WallData>("MassSpringWallData");
    AuxinGradient::CellDataAttr &cellADAttr = mesh->attributes().attrMap<CCIndex,AuxinGradient::CellData>("AuxinCellData");

    for(CCIndex e : cs.edges()) {
      MassSpring::WallData &eD = wallAttr[e];
  
      auto eb = cs.edgeBounds(e);
      double length = norm(indexAttr[eb.first].pos - indexAttr[eb.second].pos);
      double strain = (length - eD.restLength)/eD.restLength;   
     
       for(CCIndex e : cs.edges()) {
        eD.wallGrowthFactor = 0;
        int numIncidFaces = 0;
        for(CCIndex face: cs.incidentCells(e,2))
        {
          ++numIncidFaces;
          //CellData &cD = auxinGridAttr[face];
          eD.wallGrowthFactor += (cellADAttr)[face].a;
        }
        eD.wallGrowthFactor *=1./numIncidFaces;
        mdxInfo << "Edge rest length before growth is: " <<   eD.restLength << endl;
  	if(strain > strainThresh){ 
               eD.restLength += eD.restLength*growthDt*strain * growthRate;
               if(eD.restLength > length) {
                mdxInfo << "Restlength greater than length, rate too high?" << endl;
                eD.restLength = length;
               }
        }

         mdxInfo << "Growth factors multiplied: " << growthDt*eD.wallGrowthFactor <<endl;
         //mdxInfo << "wall growth factor is: " <<  eD.wallGrowthFactor << endl;
         //eD.restLength +=  eD.restLength*growthDt*eD.wallGrowthFactor;   
         mdxInfo<< "Edge rest length after growth is: "<< eD.restLength << endl;
         mdxInfo << " " << endl;
  	      //eD.restLength += eD.restLength * strain * growthDt * wallGrowth* eD.wallGrowthFactor; 
          //if(eD.restLength > length) {
          //  mdxInfo << "Restlength greater than length, rate too high?" << endl;
          //  eD.restLength = length;
          //}
        //}
      }           
     // if(strain > strainThresh) {
     //   eD.restLength += eD.restLength * strain * growthDt * growthRate; 
     //   if(eD.restLength > length) {
     //     mdxInfo << "Restlength greater than length, rate too high?" << endl;
     //     eD.restLength = length;
      //  }
     // }
    }           
    return true;
  }*/


  void MassSpring::Subdivide::splitCellUpdate(Dimension dim, const CCStructure &cs, 
      const CCStructure::SplitStruct& ss, CCIndex otherP, CCIndex otherN,  double interpPos)
  {
    if(dim == 1) {
      if(!wallAttr or !vtxAttr)
        std::cout << "MassSpringSubdivide::splitCellUpdate Wall attribute not set" << std::endl;
      else {
        double restLength = (*wallAttr)[ss.parent].restLength;
        (*wallAttr)[ss.childP].restLength = interpPos * restLength;
        (*wallAttr)[ss.childN].restLength = (1.0 - interpPos) * restLength;
        (*wallAttr)[ss.childP].wallStiffnessAssigned = (*wallAttr)[ss.parent].wallStiffnessAssigned;
        (*wallAttr)[ss.childN].wallStiffnessAssigned = (*wallAttr)[ss.parent].wallStiffnessAssigned;
        (*wallAttr)[ss.childP].wallStiffness = (*wallAttr)[ss.parent].wallStiffness;
        (*wallAttr)[ss.childN].wallStiffness = (*wallAttr)[ss.parent].wallStiffness;
        (*wallAttr)[ss.childP].noGrowth = (*wallAttr)[ss.parent].noGrowth;
        (*wallAttr)[ss.childN].noGrowth = (*wallAttr)[ss.parent].noGrowth;


        auto Pbounds = cs.edgeBounds(ss.childP);
        auto Nbounds = cs.edgeBounds(ss.childN);
        //CCIndex firstVtx; 
        //CCindex secondVts;
        Point3u dirichletFirst;
        Point3u dirichletSecond;
        Point3u dirichletNew = Point3u(0, 0, 0);
        
        if (Pbounds.first == ss.membrane)
             dirichletFirst = (*vtxAttr)[Pbounds.second].dirichlet;
        else 
             dirichletFirst = (*vtxAttr)[Pbounds.first].dirichlet;
        if (Nbounds.first == ss.membrane)
             dirichletSecond = (*vtxAttr)[Nbounds.second].dirichlet;
        else 
             dirichletSecond = (*vtxAttr)[Nbounds.first].dirichlet;

        if (dirichletFirst.x() && dirichletSecond.x())
            dirichletNew.x() = 1;
        if (dirichletFirst.y() && dirichletSecond.y())
            dirichletNew.y() = 1;
        if (dirichletFirst.z() && dirichletSecond.z())
            dirichletNew.z() = 1;
        (*vtxAttr)[ss.membrane].dirichlet = dirichletNew;
        
      }
    }
    if(dim == 2) {
      if(!cellModelAttr)
        std::cout << "MassSpring::Subdivide::splitCellUpdate Cell attribute not set" << std::endl;
      else {
        MassSpring::CellModelData &pCD = (*cellModelAttr)[ss.childP];
        MassSpring::CellModelData &nCD = (*cellModelAttr)[ss.childN];
        MassSpring::CellModelData &rCD = (*cellModelAttr)[ss.parent];

        pCD = nCD = rCD; // just copy all fields
        // What to do with PIN? For gradient it gets recalculated at every step.
      }
     }		
  }

  //// Render for polarizer vector field from gradient of the  diffusion of a signal
  bool PolarizerRender::initialize(QWidget *parent)
  {
    // Get the tissue and cell graph names
    SourceCC = currentMesh()->ccName();
    if(SourceCC.isEmpty())
      throw(QString("No input cell complex specified"));
    OutputCC = parm("Output CC");
    DrawPolarizer = stringToBool(parm("Draw Polarizer"));
    PolarizerSize = parm("Polarizer Vector Size").toDouble();
    
    
    mdxInfo << "The polarizer render has been initialized" << endl;
    return true;
  }
  
  bool PolarizerRender::run(Mesh *mesh)
  {
    // Get the cell complexes
    if(!mesh->exists(SourceCC))
      throw(QString("Specified input cell complex does not exist"));
    CCStructure &src = mesh->ccStructure(SourceCC);
    CCStructure &out = mesh->ccStructure(OutputCC);
    cellDAttr = &mesh->attributes().attrMap<CCIndex,MassSpring::CellModelData>("CellModelData");

    out = CCStructure(2);

    // We also need the attribute map
    CCIndexDataAttr &indexAttr = mesh->indexAttr();
    mdxInfo << " I could access the mesh attributes" << endl;
    if(DrawPolarizer) { 
      // Put the neighborhoods into the cell complex
      forall(CCIndex f, src.cellsOfDimension(2)) {
        mdxInfo << " I can make the faces loops" << endl;
        CCIndexData &V = indexAttr[f];
        mdxInfo << " I got the faces attributes as CCIndex " << endl;
        MassSpring::CellModelData &cCD = (*cellDAttr)[f];
        mdxInfo << " I could access the sprinModel cell attributes" << endl;
        if(DrawPolarizer) {
          // Add the vertices
          CCIndex v1 = CCIndexFactory.getIndex();
          out.addCell(v1);
          CCIndex v2 = CCIndexFactory.getIndex();
          out.addCell(v2);

          CCIndexData &V1 = indexAttr[v1];
          CCIndexData &V2 = indexAttr[v2];

          V1.pos = V.pos;
          V2.pos = V.pos + cCD.polarizerDir * PolarizerSize;

          // Add the edge
          CCIndex e = CCIndexFactory.getIndex();
          out.addCell(e, +v1 -v2);
        }
       }
    }
    mesh->drawParms(OutputCC).setGroupVisible("Vertices", true);
    mesh->drawParms(OutputCC).setGroupVisible("Edges", true);
    //mesh->drawParms(OutputCC).setGroupVisible("Faces", true);

    mesh->updateAll(OutputCC);

    return true;
  };

  bool VisualizeAssignedCellProperties::run(Mesh &mesh, const MassSpring::CellModelAttr &cellSpringAttr, const MassSpring::VertexModelAttr &vertexSpringAttr, const MassSpring::WallDataAttr &wallSpringAttr, const QString &visualize, const QString &signal)
  {
      if(cellSpringAttr.size() == 0) {
        mdxInfo << QString("%1::run Cell Spring attribute is empty").arg(name()) << endl;
        return true;
      }
      //if(wallSpringAttr.size() == 0) {
      //  mdxInfo << QString("%1::run Wall Spring attribute is empty").arg(name()) << endl;
      //  return true;
      //}

      // Get the signal attribute
      auto &signalAttr = mesh.signalAttr<double>(signal);
      //CCSignedIndexDoubleAttr *signalEdgeAttr = 0; // = currentMesh()->attributes().attrMap<CCSignedIndex, double>(signal)

      //*signalAttr = 0;
      // The first time it is created we'll set the color map and bounds
      //if(mesh.signalColorMap(signal).size() == 0)
       mesh.signalColorMap(signal).setColors("Jet");


      if(visualize == "Cell Type") {
        for(auto pr : cellSpringAttr) {
          signalAttr[pr.first] = pr.second.cellTypeSpr;
        }
        mesh.setSignalUnit("Cell Type", signal);
      } else if(visualize == "Max Area2Division") {
        for(auto pr : cellSpringAttr) {
          signalAttr[pr.first] = pr.second.maxAreaDiv;
        }
      mesh.setSignalUnit("Max Area2Division", signal);
      }
      else if(visualize == "Cell Parall Growth") {
        for(auto pr : cellSpringAttr) {
          signalAttr[pr.first] = pr.second.cellParallGrowth;
        }
      mesh.setSignalUnit("Cell Parall Growth", signal);
      }
      else if(visualize == "Cell Ortho Growth") {
        for(auto pr : cellSpringAttr) {
          signalAttr[pr.first] = pr.second.cellOrthoGrowth;
        }
      mesh.setSignalUnit("Cell Ortho Growth", signal);
      }
      else if(visualize == "Cell Iso Growth") {
        for(auto pr : cellSpringAttr) {
          signalAttr[pr.first] = pr.second.cellIsoGrowth;
        }
      mesh.setSignalUnit("Cell Iso Growth", signal);
      }
      else if(visualize == "Cell Pressure") {
        for(auto pr : cellSpringAttr) {
          signalAttr[pr.first] = pr.second.cellPressure;
        }
      mesh.setSignalUnit("Cell Pressure", signal);
      }
      else if(visualize == "Strain Based Growth Rate") {
        for(auto pr : cellSpringAttr) {
          signalAttr[pr.first] = pr.second.strBasedGrRate;
        }
      mesh.setSignalUnit("Strain Based Growth Rate", signal);
      }
      else if(visualize == "Growth Rate Strain Parall") {
        for(auto pr : cellSpringAttr) {
          signalAttr[pr.first] = pr.second.growthRateStrParall;
        }
      mesh.setSignalUnit("Growth Rate Strain Parall", signal);
      }
      else if(visualize == "Growth Rate Strain Ortho") {
        for(auto pr : cellSpringAttr) {
          signalAttr[pr.first] = pr.second.growthRateStrOrtho;
        }
      mesh.setSignalUnit("Growth Rate Strain Ortho", signal);
      }
      else if(visualize == "Strain treshold") {
        for(auto pr : cellSpringAttr) {
          signalAttr[pr.first] = pr.second.strainThreshold;
        }
      mesh.setSignalUnit("Strain treshold", signal);
      }
      else if(visualize == "Special Growth Concentration") {
        for(auto pr : cellSpringAttr) {
          signalAttr[pr.first] = pr.second.specialGrCon;
        }
      mesh.setSignalUnit("Special Growth Concentration", signal);
      }
      /*else if(visualize == "Edge Stiffness")
      {
        signalEdgeAttr = &currentMesh()->attributes().attrMap<CCSignedIndex, double>(signal);
        for(auto pr : wallSpringAttr) {
          (*signalEdgeAttr)[pr.first] = pr.second.wallStiffness;
        }
      mesh.setSignalUnit("Edge Stiffness", signal);


      }*/
      else
        throw QString("%1::run Invalid visualize choice %2").arg(name()).arg(visualize);

      return true;

  }
  /*
   * Cell Disk
   */

  // Initialize the main solver process
  bool CellDisk::initialize(QWidget *parent)
  {
    mesh = currentMesh();
    if(!mesh)
      throw(QString("CellDiskAuxin::initialize No current mesh"));

    // Initialize the tissue
    if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw(QString("CellDisk::initialize Cannot make Tissue Process:" + parm("Tissue Process")));
    tissueProcess->initialize(parent);;

    // Initialize the solver derivatives provider
    if(!getProcess(parm("Growth Signal Process"), auxinProcess))
      throw(QString("CellDisk::initialize Unable to make Growth Signal Process:") + parm("Growth Signal Process"));
    auxinProcess->initialize(parent);
   
    // Initialize the Get wall stiffness process
    if(!getProcess(parm("Set Stiffness Gradient"), setWallStiffness))
      throw(QString("CellDisk::initialize Unable to set wall stiffness gradient:") + parm("Set Stiffness Gradient"));
    setWallStiffness->initialize(parent);
     
    
 
    // Initialize the solver derivatives provider for polarizer assigment
    if(!getProcess(parm("Polarizer Process"), polarizerProcess))
      throw(QString("CellDisk::initialize Unable to make Polarizer Process:") + parm("{Polarizer Process"));
    polarizerProcess->initialize(parent);

    if(!getProcess(parm("Mass Spring Process"), massSpringProcess))
      throw(QString("CellDisk::initialize Cannot make Mass Spring Process:") + parm("Mass Spring Process"));
    massSpringProcess->initialize(parent); 
    
    if(!getProcess(parm("Divide Process"), divideProcess))
      throw(QString("CellDisk::initialize Cannot make Divide Process:") + parm("Divide Process"));
    divideProcess->initialize(parent);

    if(!getProcess(parm("Growth Process"), growthProcess))
      throw(QString("CellDisk::initialize Cannot make Growth Process:") + parm("Growth Process"));
    growthProcess->initialize(parent);

    // Get the split edges process
    if(!getProcess(parm("Split Edges Process"), splitEdgesProcess))
      throw(QString("CellDisk::initialize Cannot make Split Edges Process:") + parm("Split Edges Process"));	 

    // Get the split edges process
    if(!getProcess(parm("Only L1 Periclinal Growht Process"), onlyPeriL1grProcess))
      throw(QString("CellDisk::initialize Cannot make Only L1 Periclinal Growht Process:") + parm("Only L1 Periclinal Growht Process"));	
    maxAuxIter = parm("Max Auxin Iterations").toDouble();
    maxPolIter = parm("Max Polarizer Iterations").toDouble();
    maxMassSpringIter = parm("Max Mass Spring Iterations").toDouble();


    savingTimes = parm("Simulation time to save mesh");
    listSavingTimes = savingTimes.split(',');
    
    if(!getProcess(parm("Save Mesh Process Name"), meshSave))
        throw(QString("FemMembranes::initialize Unable to make save mesh process : %1").arg(parm("Save Mesh Process Name")));
    //meshSave->initialize(parent);
    selfStopTime = parm("Self Stop Time").toInt();
    //get the first polarizer simulation done
    //int solverIter = 0;

    //bool polarizerDifNotConvergedInit = true;
    //while(polarizerDifNotConvergedInit and solverIter < maxPolIter){
    //    polarizerDifNotConvergedInit = polarizerProcess->step();
    //    solverIter++;
    //} 
    return true;
  }

  bool CellDisk::rewind(QWidget *parent)
  {
    // To rewind, we'll reload the mesh
    mesh = currentMesh();
    if(!mesh or mesh->file().isEmpty())
      throw(QString("CellDisk::rewind No current mesh, cannot rewind"));
    MeshLoad meshLoad(*this);
    meshLoad.setParm("File Name", mesh->file());
    return meshLoad.run();
  }

  bool CellDisk::step()
  {
    //static int stepCount = 0;
    int solverIter = 0;
    bool massSpringNotConverged = true;
    bool polarizerDifNotConverged = true;
    bool auxinNotConverged = true; 
    // Do auxin simulation
    if(!auxinProcess)
      throw(QString("CellDisk::step Growth Signal Process is null"));
    // if auxin is done, get the polarizertropic diffusion process, but only if step is a multiple of 10
    while(auxinNotConverged and solverIter<maxAuxIter){
      auxinNotConverged = auxinProcess->step();
      solverIter++; 
    }
    if(!auxinNotConverged /*and stepCount%3== 0*/)
    {
      //assign the differential stiffness if you have to
      if(!setWallStiffness)   
         throw(QString("CellDisk::step set wall stiffness Process is null"));
      setWallStiffness->step();
      // Do polarizertropic diffusion simulation
      if(!polarizerProcess)
       throw(QString("CellDisk::step Polarizer Process is null"));
      // if polarizertropic process is done, go to the mechanics of mass spring
      solverIter = 0;
      while(polarizerDifNotConverged and solverIter < maxPolIter){
        polarizerDifNotConverged = polarizerProcess->step();
        solverIter++;  
      }
      if(!polarizerDifNotConverged)
      {
       mdxInfo<< "Equilibrium reached for Polarizer" << endl;
       mdxInfo << " " << endl;
       if(!massSpringProcess)
        throw(QString("CellDisk::step Mass Spring Process is null"));
        solverIter = 0;
        while(massSpringNotConverged and solverIter < maxMassSpringIter) {
          massSpringNotConverged = massSpringProcess->step();
          solverIter++;
        }
        //and now do the rest, which is not requiring any convergence check
        if (!massSpringNotConverged)
        {
         
          static int screenShotCount = 0;
          // Take screenshot
          QString fileName = QString("OvuleGrowth-%1.JPG").arg(screenShotCount++, 4, 10, QChar('0'));
          takeSnapshot(fileName);
          for (int i=0; i<listSavingTimes.size(); ++i)
          {
            if(listSavingTimes[i].toInt() == screenShotCount)
            {
              meshName = parm("Name of mesh");
              QString fileNameMesh = QString("%1-%2.mdxm").arg(meshName).arg(screenShotCount, 4, 10, QChar('0'));
              meshSave->run(mesh, fileNameMesh, true);
              break;
            }
          }
          mdxInfo << "Equilibrium reached for the mass spring mechanics" << endl;
          mdxInfo << " " << endl;
          if ( screenShotCount > selfStopTime)
          {
            mdxInfo << "I have reached the stopping condition" << endl;
            throw(QString("Reached maximal growth iteration specified by user"));
            return true;
 
          }
          //propagate the noL1Growth to membrane edges of divided cells
          if(!onlyPeriL1grProcess)
            throw(QString("CellDisk::step Only L1 Periclinal growht Process is null"));
          onlyPeriL1grProcess->run();
          if(!growthProcess)
           throw(QString("CellDisk::step Growth Process is null"));
          growthProcess->step();

          // Split large-enough cells
          if(!divideProcess)
           throw(QString("CellDisk::step Divide Process is null"));
          divideProcess->step();
          auxinProcess->initSolver(&tissueProcess->tissue().dualGraph());
          massSpringProcess->initSolver(&tissueProcess->tissue().cellStructure());
          polarizerProcess->initSolver(&tissueProcess->tissue().dualGraph());

          // Split edges
          //if(!splitEdgesProcess)
          // throw(QString("CellDisk::step Split Edges Process is null"));
          //splitEdgesProcess->run(*mesh, tissueProcess->tissue().cellStructure(), 
          //                    splitEdgesProcess->parm("Max Length").toDouble(), &divideProcess->subdivider());
        }
      }
    }
    //do mass spring anyway also if you are not in the multiple of 10 step
    /*else if(!auxinNotConverged)
    {
     // Do mass spring and if not converged, return
     if(!massSpringProcess)
      throw(QString("CellDisk::step Mass Spring Process is null"));
     solverIter = 0;
     while(massSpringNotConverged and solverIter < maxMassSpringIter){
        massSpringNotConverged = massSpringProcess->step();
        solverIter++;
     } 
     //if the mass spring has converged, move on with the rest
     if (!massSpringNotConverged){
      //save a snapshot of the simulation
      static int screenShotCount = 0;
      // Take screenshot
      QString fileName = QString("OvuleGrowth-%1.JPG").arg(screenShotCount++, 4, 10, QChar('0'));
      takeSnapshot(fileName);

      mdxInfo << "Equilibrium reached for the mass spring mechanics" << endl;
      mdxInfo << " " << endl;
      // Grow
      if(!growthProcess)
       throw(QString("CellDisk::step Growth Process is null"));
      growthProcess->step();

      //Split large-enough cells
      if(!divideProcess)
       throw(QString("CellDisk::step Divide Process is null"));
      divideProcess->step();
   
      auxinProcess->initSolver(&tissueProcess->tissue().dualGraph());
      polarizerProcess->initSolver(&tissueProcess->tissue().dualGraph());
      massSpringProcess->initSolver(&tissueProcess->tissue().cellStructure());

      //Split edges
      //if(!splitEdgesProcess)
      // throw(QString("CellDisk::step Split Edges Process is null"));
      //splitEdgesProcess->run(*mesh, tissueProcess->tissue().cellStructure(), 
      //                        splitEdgesProcess->parm("Max Length").toDouble(), &divideProcess->subdivider());
     }
    }*/
    //stepCount++;

    // Need a way to determine if any subdivision above actually divided anything, for now assume it did
    //auxinProcess->initSolver(&tissueProcess->tissue().dualGraph());
    //polarizerProcess->initSolver(&tissueProcess->tissue().dualGraph());
    //massSpringProcess->initSolver(&tissueProcess->tissue().cellStructure());
 
    return true;
  }

  // Initialize to grab subdivider
  bool CellDiskDivide::initialize(QWidget *parent)
  {
    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("CellDiskDivide::step No current mesh"));

    // Call base initialize
    if(!CellTissueCellDivideOvule::initialize(parent))
      return false;

    if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw(QString("CellDisk::initialize Cannot make tissue process"));
    tissueProcess->initialize(parent);


    //propagationDirection = stringToPoint3d(parm("Preferred Direction propagation"));
    
    // Setup subdivision objects
    subdiv.mdx = MDXSubdivide(*mesh);
    subdiv.massSpring = MassSpring::Subdivide(mesh->attributes().attrMap<CCIndex,MassSpring::WallData>("MassSpringWallData"), mesh->attributes().attrMap<CCIndex,MassSpring::CellModelData>("CellModelData"),  mesh->attributes().attrMap<CCIndex,MassSpring::VertexModelData>("VertexModelData"));
    subdiv.auxinGradient = AuxinGradient::Subdivide(tissueProcess->tissue().indexAttr(), mesh->attributes().attrMap<CCIndex,AuxinGradient::CellData>("AuxinCellData"), 
                                                    mesh->attributes().attrMap<CCSignedIndex,AuxinGradient::EdgeData>("AuxinEdgeData") 
                                                    /*, mesh->attributes().attrMap<CCIndex,MassSpring::CellModelData>("CellModelData")*/);
    return true;
  }

  void CellDiskSubdivide::splitCellUpdate(Dimension dim, const CCStructure &cs, const CCStructure::SplitStruct& ss, 
          CCIndex otherP, CCIndex otherN, double interpPos)
  {
    mdx.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
    massSpring.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
    auxinGradient.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
  }

  // Run a step of cell division
  bool CellDiskDivide::step() 
  { 
    // Pass our subdivide
    return CellTissueCellDivideOvule::step(currentMesh(), &subdiv);
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //copy here the modified MDXProcessCellDivide
   // Divide this so it can be called on any tissue
  bool CellTissueCellDivideOvule::initialize(QWidget* parent)
  {
    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("CellTissueCellDivideOvile::initialize No current mesh"));

    processParms();

    // Initialize the tissue
    if(!getProcess(TissueParmsProcessName, cellTissueProcess))
      throw(QString("CellTissueCellDivideOvule::initialize: Unable to create CellTissueProcess: %1").arg(TissueParmsProcessName));
    cellTissueProcess->initialize(*mesh);
    tissue = &cellTissueProcess->tissue();
    if(!tissue)
      throw(QString("CellTissueCellDivideOvule::initialize: Unable to get cell tissue"));

    TissueName = cellTissueProcess->parm("Tissue");
    TissueDualName = cellTissueProcess->parm("Tissue Dual");

    // Set the cell division parms
    //tissue->processCellDivideParms(*this);

    // Just use the basic MDX data for subdivision
    subdiv = MDXSubdivide(*mesh);

    return true;
  }

  bool CellTissueCellDivideOvule::processParms()
  {
    // Process parameters
    CellMaxArea = parm("Cell Max Area").toDouble();
    Verbose = stringToBool(parm("Verbose"));
    TissueParmsProcessName = parm("Tissue Process");
    DivAlg = parm("Division Algorithm").toInt();
    selectivePropagation =stringToBool(parm("Selective propagation"));
    return true;
  }

  // Run a step of cell division
  bool CellTissueCellDivideOvule::step(Mesh *mesh, Subdivide *subdiv)
  {
    if(!mesh)
      throw(QString("CellTissueCellDivideOvule.step: Mesh not initialized"));
    if(!tissue)
      throw(QString("CellTissueCellDivideOvule.step: Tissue not initialized"));

    CCStructure &cs = cellTissueProcess->tissue().cellStructure();

    CCIndexDataAttr &indexAttr = tissue->indexAttr();
    /*MassSpring::CellModelAttr &*/cellModelAttr = &mesh->attributes().attrMap<CCIndex,MassSpring::CellModelData>("CellModelData");
    /*AuxinGradient::CellDataAttr &*/cellDataAttr = &mesh->attributes().attrMap<CCIndex,AuxinGradient::CellData>("AuxinCellData");

    // Check if we need to use selection
    // implement later

    // Find cells to divide
    std::vector<CCIndex> D;

    for(CCIndex c : cs.faces()) {
      MassSpring::CellModelData &cCD = (*cellModelAttr)[c];
      AuxinGradient::CellData &cDA = (*cellDataAttr)[c];
      CCIndexData &cIdx = indexAttr[c];
      CellMaxAreaOvule = &cCD.maxAreaDiv;
      //cDA.propFixConDivision = false;
      double area = indexAttr[c].measure;
      if ( *CellMaxAreaOvule<=1.e-6 or cCD.cellTypeSpr == MassSpring::NormalCell)
        *CellMaxAreaOvule = CellMaxArea;
      if (area>*CellMaxAreaOvule){
         D.push_back(c);
           if (cDA.fixedConc == true)
           {
              mdxInfo << "I am handling the concentration propagation" << endl;
              //this is a pseudo position of nearest fixed cell: if among the neighbour cells with fixed concentration there is an L1, this wins
              cDA.posNearestFixed = Point3d(1.e30, 1.e30, 1.e30);
              double distance = 1.e20; 
              int numFixed = 0;
              for(CCIndex cf : cs.neighbors(c))
              {  
                 AuxinGradient::CellData &cfDA = (*cellDataAttr)[cf];
                 CCIndexData &cfIdx = indexAttr[cf];
                 MassSpring::CellModelData &cmA = (*cellModelAttr)[cf];
                 if (cfDA.fixedConc == true){
                   if ((*cellModelAttr)[c].cellTypeSpr != MassSpring::CellTypeSpr::L1)
                     numFixed++;
                   else if ((*cellModelAttr)[cf].cellTypeSpr == MassSpring::CellTypeSpr::L1)
                     numFixed++;
                   if (cmA.cellTypeSpr == MassSpring::CellTypeSpr::L1 and (*cellModelAttr)[c].cellTypeSpr != MassSpring::CellTypeSpr::L1)
                   {
                      cDA.posNearestFixed = cfIdx.pos;
                      distance = -1.e30;
                      mdxInfo << "I am in the case of being connected to a L1 fixed source" << endl;
                      //continue;
                   }
                   double distanceTemp = (cIdx.pos - cfIdx.pos).norm();
                   if (distanceTemp < distance)
                   {
                      distance = distanceTemp;
                      cDA.posNearestFixed = cfIdx.pos;
                   }
                 }
              }  
              if (numFixed > 1 or !selectivePropagation ) //n.b for epidermal growth models this number should be 2!!!!!
                cDA.propFixConDivision = true;
              
              
           }
         cDA.copyPolDir = cCD.polarizerDir;
      }
    }

    if(D.empty())
      return false;
    
    
    forall(const CCIndex &cell, D) {
      if(Verbose)
        mdxInfo << "Dividing cell (" << to_string(cell).c_str() << "), with area " << indexAttr[cell].measure << endl;
      if ((*cellModelAttr)[cell].divAlg == 0 )
       //tissue->divideCell(cell, subdiv);
       divideCell2d(tissue->cellStructure(), indexAttr, cell, *this, subdiv);
      //else{ 
        //CellTissue::DivAlg divAlgo = CellTissue::ASSIGNED_VECTOR_TRHOUGH_CENTROID;
        //tissue->divideCell(cell, subdiv, divAlgo, (cellAttr)[cell].divVector);
      //}
    }
    // Re create dual and update geometry
    tissue->createDual();
    tissue->updateGeometry();

    // Update mesh points, edges, surfaces
    mesh->updateAll(tissue->tissueName());
    mesh->updateAll(tissue->tissueDualName());

    return true; 
  }
  

  REGISTER_PROCESS(AuxinSolver);
  REGISTER_PROCESS(AuxinGradient);
  REGISTER_PROCESS(AuxinRender);
  REGISTER_PROCESS(AuxinSetCellType);
  REGISTER_PROCESS(AuxinResetCellType);

  REGISTER_PROCESS(PolarizerSetFixedConc);
  REGISTER_PROCESS(PolarizerResetFixedConc);
  REGISTER_PROCESS(PolarizerDiffusionSolver);
  REGISTER_PROCESS(PolarizerDiffusion);
  REGISTER_PROCESS(PolarizerRender);


  REGISTER_PROCESS(MassSpringSolver);
  REGISTER_PROCESS(MassSpring);
  REGISTER_PROCESS(MassSpringGrowth);
  REGISTER_PROCESS(SetDirichlet);
  REGISTER_PROCESS(SetWallStiffness);
  REGISTER_PROCESS(ResetWallStiffness);
  REGISTER_PROCESS(SpringSetCellType);
  REGISTER_PROCESS(SpringResetCellType);
  REGISTER_PROCESS(OnlyPericlinalGrowthL1);
  REGISTER_PROCESS(VisualizeAssignedCellProperties);
  REGISTER_PROCESS(CellDiskVisDirichlet);
  REGISTER_PROCESS(CellDisk);
  REGISTER_PROCESS(CellDiskTissue);
  REGISTER_PROCESS(CellDiskDivide);
  REGISTER_PROCESS(CellDiskSplitEdges);
}
