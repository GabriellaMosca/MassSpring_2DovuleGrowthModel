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

#ifndef CELL_DISK_HPP
#define CELL_DISK_HPP

#include <Solver.hpp>
#include <Attributes.hpp>
#include <MeshProcessSystemRender.hpp>
#include <MeshProcessStructure.hpp>
#include <MeshProcessSystem.hpp>
#include <MDXProcessTissue.hpp>
#include <MDXProcessCellDivide.hpp>
#include <MDXProcessFem.hpp>
#include <MeshUtils.hpp>
#include <Fem.hpp>
using namespace mdx;
using namespace fem;
namespace CellDisk
{
  /*
   * Auxin simulation, up the gradient
   */
  class AuxinGradient : public Process, public SolverDerivs<Point1d>
  {
  public:
    // Define cell type
    enum CellType {GroundCell, SourceCell, SinkCell, NumCellTypes};
    double maxA, minA; 
    static CellType stringToCellType(const QString &str)
    {
      if(str == "Ground Cell")
        return(GroundCell);
      else if(str == "Source Cell")
        return(SourceCell);
      else if(str == "Sink Cell")
        return(SinkCell);
      else 
        throw(QString("Bad cell type %1").arg(str));
    }
  
    // Structure to store the cell data
    struct CellData
    { 
      double a = 0; // Morphogen concentration
      double perimeter = 0;
      double area = 0;
  
      std::vector<Point3d *> posVec;
      CCIndexVec vtxVec;
  
      double amount = 0; // Amount for source or sink
      uint cellType = GroundCell; // Type of cell
      bool fixedConc = false; // Fixed concentration, or prod/decay
      bool propFixConDivision = false;      
      Point3d posNearestFixed = Point3d(1.e30, 1.e30, 1.e30);
      bool propagateLaterally = false;
      Point3d copyPolDir = Point3d(0., 0., 0.); 
      CellData() {}
  
      bool operator==(const CellData &other) const
      {
        if(a == other.a and perimeter == other.perimeter
            and area == other.area and amount == other.amount 
            and cellType == other.cellType and fixedConc == other.fixedConc and propFixConDivision ==other.propFixConDivision and posNearestFixed == other.posNearestFixed and propagateLaterally == other.propagateLaterally and copyPolDir == other.copyPolDir)
          return true;
        return false;
      }
    };
    typedef AttrMap<CCIndex,CellData> CellDataAttr;
    
    // Structure to store the edge data
    struct EdgeData
    {
      double pin; // Morphogen concentration
      double flux;
  
      EdgeData() : pin(0), flux(0) {}
  
      bool operator==(const EdgeData &other) const
      {
        if(pin == other.pin and flux == other.flux)
          return true;
        return false;
      }
    };
    typedef AttrMap<CCSignedIndex,EdgeData> EdgeDataAttr; 

    class Subdivide : virtual public mdx::Subdivide
    {
    public:
      Subdivide() {}
      Subdivide(CCIndexDataAttr &_indexAttr, CellDataAttr &_cellAttr,  EdgeDataAttr &_edgeAttr)
      {
        cellAttr = &_cellAttr;
        edgeAttr = &_edgeAttr;
        indexAttr = &_indexAttr;
        //propagVector = &_propagVector;
        //cellModelAttr = &_cellModelAttr;
      }

      void splitCellUpdate(Dimension dim, const CCStructure &cs, const CCStructure::SplitStruct& ss, 
          CCIndex otherP = CCIndex(), CCIndex otherN = CCIndex(),  double interpPos = 0.5);

      CellDataAttr *cellAttr = 0;
      EdgeDataAttr *edgeAttr = 0;
      CCIndexDataAttr  *indexAttr = 0;
      Point3d *propagVector = 0;
      //MassSpring::CellModelAttr *cellModelAttr = 0;

    };

    // Constructor
    AuxinGradient(const Process &process) : Process(process) 
    {
      setName("Model/Cell Ovule Growth/10 Growth Signal/b Growth Signal Gradient");
      setDesc("Transport on a cellular grid");
      setIcon(QIcon(":/images/CellDisk.png"));

      insertParm("Production", "Production coefficient for ground cells", "10.0", 0);
      insertParm("Transport", "Transport coefficient", "500.0", 1);
      insertParm("Diffusion", "Diffusion coefficient", "5.0", 2);
      insertParm("Decay", "Decay coefficient", "0.1", 3);
    }
    bool initialize(QWidget *parent) { return true; } // Should not be run from the GUI
    bool initialize(CellTissue &tissue, CCIndexDataAttr &indexAttr, CellDataAttr &cellAttr, EdgeDataAttr &edgeAttr);

    // Reimplemented solver methods, Jacobian is tricky for up-the-gradient PIN
    void setValues(const SolverT &solver, CCIndex c, const VectorT &values)
    {
      CellData &cCD = (*cellAttr)[c];
      cCD.a = values[0];
    }
    void getValues(const SolverT &solver, CCIndex c, VectorT &values)
    {
      CellData &cCD = (*cellAttr)[c];
      values[0] = cCD.a;
    }    
    void calcDerivatives(const SolverT &solver, VertexAttr &vertexData);
    void initDerivs(SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr);

  private:
    CellTissue *tissue = 0;
    CCIndexDataAttr *indexAttr = 0;
    CellDataAttr *cellAttr = 0;
    EdgeDataAttr *edgeAttr = 0;

    // Parameters
    double Prod = 0;
    double Trans = 0;
    double Diff = 0;
    double Decay = 0;
  };
 
  // Read/write Cell and Edge data
  bool inline readAttr(AuxinGradient::CellData &m, const QByteArray &ba, size_t &pos) 
  { 
    return mdx::readChar((char *)&m, sizeof(AuxinGradient::CellData), ba, pos);
  }
  bool inline writeAttr(const AuxinGradient::CellData &m, QByteArray &ba)
  { 
    return mdx::writeChar((char *)&m, sizeof(AuxinGradient::CellData), ba);
  }
  bool inline readAttr(AuxinGradient::EdgeData &m, const QByteArray &ba, size_t &pos) 
  { 
    return mdx::readChar((char *)&m, sizeof(AuxinGradient::EdgeData), ba, pos);
  }
  bool inline writeAttr(const AuxinGradient::EdgeData &m, QByteArray &ba)
  { 
    return mdx::writeChar((char *)&m, sizeof(AuxinGradient::EdgeData), ba);
  }
  
  // Process to set sources and sinks
  class AuxinSetCellType : public Process
  {
  public:
    AuxinSetCellType(const Process &process) : Process(process) 
    {
      setName("Model/Cell Ovule Growth/10 Growth Signal/a Set Cell Type/Assign");
      setDesc("Set the cell type.");
      setIcon(QIcon(":/images/CellType.png"));

      addParm("Cell Type", "Identity of the cell", "Ground Cell", QStringList() << "Ground Cell" << "Source Cell" << "Sink Cell");
      addParm("Amount", "Amount of source sink", "10.0");
      addParm("Fixed Conc", "Amount is fixed concentration, or production/decay", "True", booleanChoice());
      addParm("Propagate laterally", "Propagate fixed concentration in divisions occurring with a wall orthogonal to polarization field", "False" , QStringList() << "False" << "True");
    }
    bool run();

  private:
    uint cellType;
    double Amount;
    bool FixedConc;
    bool PropagateLaterally;

  };
 
  class AuxinResetCellType : public Process
  {
  public:
    AuxinResetCellType(const Process &process) : Process(process) 
    {
      setName("Model/Cell Ovule Growth/10 Growth Signal/a Set Cell Type/Reset");
      setDesc("Set the cell type.");
      setIcon(QIcon(":/images/CellType.png"));

    }
    bool run();

  private:
    uint cellType;
    double Amount;
    bool FixedConc;
  };


  // Auxin solver
  class AuxinSolver : public Process, public Solver<Point1d>
  {
  public:
    AuxinSolver(const Process &process) : Process(process) 
    {
      setName("Model/Cell Ovule Growth/10 Growth Signal/c Growth Signal Solver");
      setDesc("Diffusion on a cellular grid");
      setIcon(QIcon(":/images/CellDisk.png"));

      addParm("Converge Threshold", "Threshold for convergence", ".1");
      addParm("Growth Signal Derivatives", "Name of process for growth signal model derivatives", "Model/Cell Ovule Growth/10 Growth Signal/b Growth Signal Gradient");
      addParm("Tissue Process", "Name of process for Cell Tissue simulation", "Model/Cell Ovule Growth/40 Cell Tissue/a Cell Tissue Process");
    }
    bool initialize(QWidget *parent);
    bool step();

    CellTissue &tissue();
    double convergeThresh = 1.e-6;
  private:
    Mesh *mesh = 0;
    
    CCIndexDataAttr *indexAttr = 0;
    AuxinGradient::CellDataAttr *cellAttr = 0;
    AuxinGradient::EdgeDataAttr *edgeAttr = 0;
    AuxinGradient *auxinProcess = 0;
    CellTissueProcess *tissueProcess = 0;
  };

  // Add a render polarity process in the GUI to make it easy to set
  class AuxinRender : public RenderPolarity
  {
  public:
    AuxinRender(const Process &process) : RenderPolarity(process)
    {
      setName("Model/Cell Ovule Growth/10 Growth Signal/d Growth Signal Render");
    }
  };

  /*
   * Mass-spring processes
   */
  class MassSpring : public Process, public SolverDerivs<Point3d>
  {
  public:
    enum CellTypeSpr {MegasporeMother, L1, SpecialCell, NormalCell, NumCellTypes };
  
    static CellTypeSpr stringToCellTypeSpr(const QString &str)
    {
      if(str == "Megaspore Mother Cell")
        return(MegasporeMother);
      else if(str == "L1 Cell")
        return(L1);
      else if(str == "Special Cell")
        return(SpecialCell);
      else if(str == "Normal Cell")
        return(NormalCell);
      else 
        throw(QString("Bad cell type %1").arg(str));
    }

    /*enum CellTypeDiffSpr {FixedConcH, FicedConcL, NumCellTypes };
  
    static CellTypeDiffSpr stringToCellTypeDiffSpr(const QString &str)
    {
      if(str == "Fixed Concentration Low")
        return(FixedConcH);
      else if(str == "Fixed Concentration High")
        return(FixedConcH);
      else 
        throw(QString("Bad cell diff spring type %1").arg(str));
    }*/

    struct CellModelData
    {
      double area = 0;
      double maxAreaDiv = 0;
      uint cellTypeSpr = NormalCell;
      bool fixDiffAnisoConc = false;
      double fixedConcVal = 0.;
      double concentration = 0.;
      //CellTissue::DivAlg divAlg = SHORTEST_WALL_THROUGH_CENTROID;
      uint divAlg = 0;
      Point3d divVector = Point3d(1., 0., 0.);
      Point3d polarizerDir = Point3d(0., 0., 0.);
      double cellParallGrowth = 1.;
      double cellOrthoGrowth = 0.;
      bool cellIsoGrowth = true;
      double cellPressure = 0.0;
      double cellStiffness = 0.0;
      double strainThreshold = 0.0;
      double strBasedGrRate = 0.0;
      double specialGrCon = 0.0;
      double growthRateStrParall = 0;
      double growthRateStrOrtho = 0;
 
      CellModelData() {}
  
      bool operator==(const CellModelData &other) const
      {
        if(area == other.area and maxAreaDiv == other.maxAreaDiv 
           and cellTypeSpr == other.cellTypeSpr and fixDiffAnisoConc == other.fixDiffAnisoConc and fixedConcVal == other.fixedConcVal and
           concentration == other.concentration and divAlg == other.divAlg and divVector == other.divVector and polarizerDir == other.polarizerDir and cellParallGrowth == other.cellParallGrowth
           and cellOrthoGrowth == other.cellOrthoGrowth and cellIsoGrowth == other.cellIsoGrowth and cellPressure == other.cellPressure and cellStiffness == other.cellStiffness and strainThreshold == other.strainThreshold and strBasedGrRate == other.strBasedGrRate and specialGrCon == other.specialGrCon and growthRateStrParall == other.growthRateStrParall and growthRateStrOrtho == other.growthRateStrOrtho)
          return true;
        return false;
      }
    };
    typedef AttrMap<CCIndex, CellModelData> CellModelAttr;


    // structure to store model vtx data
    struct VertexModelData
    {
     Point3u dirichlet = Point3u(0, 0, 0); // fixed Dirichlet boundary in x, y, z
     Point3d Neumann = Point3d(0., 0., 0.); // Set Neumann on nodes
     VertexModelData() {}
  
     bool operator==(const VertexModelData &other) const
     {
      if(dirichlet == other.dirichlet and Neumann == other.Neumann)
        return true;
      return false;
     }      
    };
    typedef AttrMap<CCIndex,VertexModelData> VertexModelAttr;
    struct WallData
    {
      double restLength = 0; // target area/length
      double wallGrowthFactor = 0; //cumulative growth factor on an edge due to signal contribution mediated over the faces the edge belongs to
      double wallStrainGrowthFactor = 0.;
      double wallStrainParallelGrowthFactor = 0.;
      double wallStrainOrthoGrowthFactor = 0.;

      double wallStrainThreshold = 0.;
      double wallStiffness = 0.;
      bool wallStiffnessAssigned = false;
      bool noGrowth = false;
      WallData() {}
  
      bool operator==(const WallData &other) const
      {
        if(restLength == other.restLength and wallGrowthFactor == other.wallGrowthFactor and wallStrainGrowthFactor == other.wallStrainGrowthFactor and wallStrainParallelGrowthFactor == other.wallStrainParallelGrowthFactor and wallStrainOrthoGrowthFactor == other.wallStrainOrthoGrowthFactor and wallStrainThreshold == other.wallStrainThreshold and wallStiffness == other.wallStiffness and wallStiffnessAssigned == other.wallStiffnessAssigned and noGrowth == other.noGrowth)
           return true;
         return false;
      }      
    };
    typedef AttrMap<CCIndex, WallData> WallDataAttr;		 

    class Subdivide : virtual public mdx::Subdivide
    {
    public:
      Subdivide() {}
      Subdivide(MassSpring::WallDataAttr &_wallAttr, MassSpring::CellModelAttr &_cellModelAttr, MassSpring::VertexModelAttr &_vtxAttr) : wallAttr(&_wallAttr), cellModelAttr(&_cellModelAttr), vtxAttr(&_vtxAttr) {}

      void splitCellUpdate(Dimension dim, const CCStructure &cs, const CCStructure::SplitStruct& ss, 
          CCIndex otherP = CCIndex(), CCIndex otherN = CCIndex(),  double interpPos = 0.5);

      MassSpring::WallDataAttr *wallAttr = 0;
      MassSpring::CellModelAttr *cellModelAttr = 0;
      MassSpring::VertexModelAttr *vtxAttr = 0;

    };

    MassSpring(const Process &process): Process(process) 
    {
      setName("Model/Cell Ovule Growth/20 Mass Spring Mechanics/e Mass Spring Parameters");

      addParm("Spring Constant", "Stiffness of the springs", "1.0");
      addParm("Turgor Pressure", "Value of the turgor pressure in the cells", "0.1");
      addParm("Tissue Process", "Name of the process for the Cell Tissue", "Model/Cell Ovule Growth/40 Cell Tissue/a Cell Tissue Process");
    }
    
    // Reimplemented SolverDerivs methods
    void setValues(const SolverT &solver, CCIndex v, const VectorT &values)
    {
      if(!indexAttr)
        throw(QString("MassSpring::setValues Attribute map not set"));
      CCIndexData &vIdx = (*indexAttr)[v];
      vIdx.pos = values;
    }
    void getValues(const SolverT &solver, CCIndex v, VectorT &values)
    {
      if(!indexAttr)
        throw(QString("MassSpring::getValues Attribute map not set"));
      const CCIndexData &vIdx = (*indexAttr)[v];
      values = vIdx.pos;
    }
    void calcDerivatives(const SolverT &solver, CCIndex v, VectorT &dx);
    void initDerivs(SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr);

    //void calcDerivativesDiff(const SolverT &solver, CCIndex v, VectorT &dx);
    //void initDerivsDiff(SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr);


    // Initialize
    bool initialize(QWidget *parent);

  private:
    Point3d calcDerivsCell(const CCStructure &cs, const CCIndexDataAttr &indexAttr, CCIndex f, CCIndex e, CCIndex v);
    Point3d calcDerivsWall(const CCStructure &cs, const CCIndexDataAttr &indexAttr, CCIndex e, CCIndex v);

    double delta = 0.;
    double springK = 0.;
    double pressure = 0.;
    //double DiffAniso = 0;


    Mesh *mesh = 0;
    CellTissueProcess *tissueProcess = 0;
    CCIndexDataAttr *indexAttr = 0;
    WallDataAttr *wallAttr = 0;
    CellModelAttr *cellDAttr = 0; 
    VertexModelAttr *vtxDAttr = 0;
  };
  
  bool inline readAttr(MassSpring::CellModelData &m, const QByteArray &ba, size_t &pos)
  {
     return mdx::readChar((char *)&m, sizeof(MassSpring::CellModelData), ba, pos);
  }
  bool inline writeAttr(const MassSpring::CellModelData &m, QByteArray &ba) 
  {
    return mdx::writeChar((char *)&m, sizeof(MassSpring::CellModelData), ba);  
  }

  bool inline readAttr(MassSpring::VertexModelData &m, const QByteArray &ba, size_t &pos)
  {
     return mdx::readChar((char *)&m, sizeof(MassSpring::VertexModelData), ba, pos);
  }
  bool inline writeAttr(const MassSpring::VertexModelData &m, QByteArray &ba) 
  {
    return mdx::writeChar((char *)&m, sizeof(MassSpring::VertexModelData), ba);
  }

  bool inline readAttr(MassSpring::WallData &m, const QByteArray &ba, size_t &pos)
  {
    return mdx::readChar((char *)&m, sizeof(MassSpring::WallData), ba, pos);
  }
  bool inline writeAttr(const MassSpring::WallData &m, QByteArray &ba) 
  {
    return mdx::writeChar((char *)&m, sizeof(MassSpring::WallData), ba);
  }


  class SetDirichlet : public Process
  {
    public:
     SetDirichlet(const Process &process) : Process(process) 
     {

      setName("Model/Cell Ovule Growth/20 Mass Spring Mechanics/a Set Dirichlet");
      setDesc("Assign Dirichlet to vertexes.");
      setIcon(QIcon(":/images/CellType.png"));

      addParm("Dirichlet", "Set Dirichlet condition on selected vtx in x,y,z", "0 0 0");
      addParm("Nodal Neumann", "Set Neumann condition on selected vtx in x,y,z", "0.0 0.0 0.0");
      addParm("Dirichlet Attribute", "Dirichlet Attribute for visualization", "Dirichlet Attribute");
     }
     bool initialize(QWidget *parent);
     bool run();
    
    private:
     Point3u dirichletFixed;
     Point3d nodalNeumann;
     fem::DirichletAttr *dirichletAttr = 0;
  };

  class SpringSetCellType : public Process
  {
     public:
     SpringSetCellType(const Process &process) : Process(process) 
     {
      setName("Model/Cell Ovule Growth/20 Mass Spring Mechanics/b Set Cell Type Mechanics and Growth/Assign");
      setDesc("Set the cell type for Mechanical simulation with cell-specific growth/mechanical parameters (they override general settings)");
      setIcon(QIcon(":/images/CellType.png"));

      addParm("Cell Type", "Identity of the cell", "L1 Cell", QStringList() <<"Megaspore Mother Cell" << "L1 Cell" << "Special Cell");
      addParm("maxArea", "Maximal area before division", "10.0");
      addParm("Isotropic Growth", "boolean for isotropic growth", "True", QStringList() <<"True" << "False" );
      addParm("Parallel Growth to Field", "Amount of growth parallel to polarizer, ranging from 0 to 1", "1.0");
      addParm("Orthogonal Growth to Field", "Amount of growth orthogonal to polarizer, ranging from 0 to 1", "1.0");
      addParm("Cell Pressure", "Specific cell pressure", "0.1");
      //addParm("Cell Stiffness", "Average cell stiffness", "100.0");
      addParm("Strain threshold", "Strain threshold for strain based growth", "0.0");
      addParm("Strain based growth rate", "Strain based growth rate", "0.0");
      addParm("Strain based growth rate par to field", "Strain based growth rate parallel to field", "0.0");
      addParm("Strain based growth rate ortho to field", "Strain based growth rate orthogonal to field", "0.0");
      addParm("Special growth rate intensity", "Growth rate specified (ignore signal diffusion), inactive if equal to 0", "0.0");
      //addParm("Division Rule", "Rule for cell division", "ORTH_MAX_GROWTH");
     }
     //bool initialize(QWidget *parent);
     bool run();

     private:
     uint cellTypeSpr;
     double maxArea;
     bool isoGrowth;
     double paralGrowth;
     double orthoGrowth;
     double cellPressure;
     //double cellStiffness; 
     double strainThresh;
     double strBasedGrRate;
     double growthRateStrParall;
     double growthRateStrOrtho;
     double specialGrCon;
  };

  class SpringResetCellType : public Process
  {
     public:
     SpringResetCellType(const Process &process) : Process(process) 
     {
      setName("Model/Cell Ovule Growth/20 Mass Spring Mechanics/b Set Cell Type Mechanics and Growth/Reset");
      setDesc("Reset the cell type for Mechanical simulation. All the parameters are assigned as to Normal Cells from the global growth and mechanics processes");
      setIcon(QIcon(":/images/CellType.png"));

     }
     //bool initialize(QWidget *parent);
     bool run();

     private:
     uint cellTypeSpr;
  };

  class SetWallStiffness : public Process
  {
    public:
     SetWallStiffness(const Process &process) : Process(process) 
     {

      setName("Model/Cell Ovule Growth/20 Mass Spring Mechanics/c Set Wall Stiffness/Assign");
      setDesc("Assign special stiffness to selected cells -->transfered then to edges.");
      setIcon(QIcon(":/images/CellType.png"));

      addParm("Wall Stiffness", "Assign stiffness to selected walls", "200.");
      addParm("Wall Stiffness based on morphogen", "Assign wall stiffness proportional to morphogen signal", "True", QStringList() << "True" << "False");
      addParm("Wall Stiffness Max", "Max wall stiffness if proportional to morphogen signal", "300" );
      addParm("Wall Stiffness Min", "Min wall stiffness if proportional to morphogen signal", "100" );
      addParm("Proportionality type", "Type of proportionality, positive or negative", "Negative", QStringList() << "Negative" << "Positive");
      addParm("Tissue Process", "Name of the process for the Cell Tissue", "Model/Cell Ovule Growth/40 Cell Tissue/a Cell Tissue Process");
      addParm("Growth Signal Gradient", "Name of Growth Signal Gradient Process", "Model/Cell Ovule Growth/10 Growth Signal/b Growth Signal Gradient");

     }
     bool initialize(QWidget *parent);
     bool run();
    
     CellTissue &tissue();
    private:
     //Mesh *mesh = 0;
 
     //CellTissueProcess *tissueProcess = 0;
     AuxinGradient *auxinGradientStiff = 0;
     double specialStiffness;
     bool stiffnessBasedMorpho;
     double maxStiffness;
     double minStiffness;
     QString proportionality; 
  };

  class ResetWallStiffness : public Process
  {
    public:
     ResetWallStiffness(const Process &process) : Process(process) 
     {

      setName("Model/Cell Ovule Growth/20 Mass Spring Mechanics/c Set Wall Stiffness/Reset");
      setDesc("Reset stiffness for all walls to the global one");
      setIcon(QIcon(":/images/CellType.png"));

      //addParm("Wall Stiffness", "Assign stiffness to selected walls", "200.");
     }
     bool run();
    
    private:
    
  };

  class OnlyPericlinalGrowthL1 : public Process
  {
    public:
     OnlyPericlinalGrowthL1(const Process &process) : Process(process) 
     {

      setName("Model/Cell Ovule Growth/20 Mass Spring Mechanics/d Only Periclinal Growth in L1");
      setDesc("Anticlinal walls in L1 are prevented from growing-- needs to be run after b Set Cell Type Mechanics and Growth/Assign");
      setIcon(QIcon(":/images/CellType.png"));

      addParm("Only Periclinal Growth in L1", "Blocks anticlinal growth for L1 later", "True", QStringList() <<"False" << "True");

     }
     //bool initialize(QWidget *parent);
     bool run();
    
     CellTissue &tissue();
    private:
     //Mesh *mesh = 0;
    bool onlyPericlinalGrL1 = false;
 
  };

  class PolarizerSetFixedConc : public Process
  {
     public:
     PolarizerSetFixedConc(const Process &process) : Process(process) 
     {
      setName("Model/Cell Ovule Growth/30 Growth Process/a Set Fixed Concentr Diffusion Polarizer/Assign");
      setDesc("Set Fixed Concentrations for Diffusive Process to set Polarizer.");
      setIcon(QIcon(":/images/CellType.png"));

      //addParm("Fixed Cell", "Set the fixed concentration-cells", "Fixed Concentration Low", QStringList() <<"Fixed Concentration High");
      addParm("Concentration", "Value of concentration the fixed concentration-cells", "1.");

     }
     //bool initialize(QWidget *parent);
     bool run();

     private:
     double fixedConcentration;
     //uint cellTypeDiffSpr;
  };

   class PolarizerResetFixedConc : public Process
  {
     public:
     PolarizerResetFixedConc(const Process &process) : Process(process) 
     {
      setName("Model/Cell Ovule Growth/30 Growth Process/a Set Fixed Concentr Diffusion Polarizer/Reset");
      setDesc("Reset Fixed Concentrations for Diffusive Process to set Polarizer.");
      setIcon(QIcon(":/images/CellType.png"));

     }
     //bool initialize(QWidget *parent);
     bool run();

     private:
     //uint cellTypeDiffSpr;
  };


  class PolarizerDiffusion : public Process, public SolverDerivs<Point1d>
  {
    public:
      PolarizerDiffusion(const Process &process): Process(process) 
      {
       setName("Model/Cell Ovule Growth/30 Growth Process/c Diffusion for Polarizer");

       addParm("Tissue Process", "Name of the process for the Cell Tissue", "Model/Cell Ovule Growth/40 Cell Tissue/a Cell Tissue Process");
       addParm("Diffusion Constant", "Value of the diffusion constant", "5.");
      }
      
    
      // Reimplemented SolverDerivs methods
      void setValues(const SolverT &solver, CCIndex c, const VectorT &values)
      {
       if(!indexAttr)
        throw(QString("PolarizerDiff::setValues Attribute map not set"));

       MassSpring::CellModelData &cCD = (*cellDAttr)[c];
       cCD.concentration = values[0];
      }
      void getValues(const SolverT &solver, CCIndex c, VectorT &values)
      {
       if(!indexAttr)
        throw(QString("PolarizerDiff::getValues Attribute map not set"));

       MassSpring::CellModelData &cCD = (*cellDAttr)[c];
       values[0] = cCD.concentration;
      }
 
      //void initDerivs(SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr);
      void calcDerivatives(const SolverT &solver, VertexAttr &vertexAttr);

      //void initDerivs(SolverT &solver, MassSpring::VertexModelAttr &vtxDAttr, MassSpring::WallDataAttr &wallAttr);

      //void calcDerivativesDiff(const SolverT &solver, CCIndex v, VectorT &dx);
      //void initDerivsDiff(SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr);


      // Initialize
      bool initialize(QWidget *parent);
      //bool initialize(CellTissue &tissue, CCIndexDataAttr &indexAttr, MassSpring::CellModelAttr &cellDAttr, MassSpring::WallDataAttr &wallAttr);

    private:
      //Point3d calcDerivsCell(const CCStructure &cs, const CCIndexDataAttr &indexAttr, CCIndex f, CCIndex e, CCIndex v);

      double delta = 1e-6;
      double DiffConst = 0.;


      Mesh *mesh = 0;
      //CellTissue *tissue = 0;
      CellTissueProcess *tissueProcess = 0;

      MassSpring::VertexModelAttr *vtxDAttr = 0;


      CCIndexDataAttr *indexAttr = 0;
      //MassSpring::WallDataAttr *wallAttr = 0;
      MassSpring::CellModelAttr *cellDAttr = 0; 
     
  };

  class PolarizerDiffusionSolver : public Process, public Solver<Point1d>
  {
    public:
      PolarizerDiffusionSolver(const Process &process): Process(process) 
      {
       setName("Model/Cell Ovule Growth/30 Growth Process/d Solver Diffusion for Polarizer");
       setDesc("Diffusion on a template for polarizer vector.");
       setIcon(QIcon(":/images/CellDisk.png"));
      

       addParm("Converge Threshold", "Threshold for convergence", ".1");
       addParm("Tissue Process", "Name of the process for the Cell Tissue", "Model/Cell Ovule Growth/40 Cell Tissue/a Cell Tissue Process");
       addParm("Polarizer Diffusion Process", "Name of Polarizer Diffusion derivatives process", "Model/Cell Ovule Growth/30 Growth Process/c Diffusion for Polarizer");
      }
      bool initialize(QWidget *parent);
      bool step();
       
      //CellTissue &tissue();

   private:
    Mesh *mesh = 0;

    CellTissueProcess *tissueProcess = 0;
    CCIndexDataAttr *indexAttr = 0;

    MassSpring::VertexModelAttr *vtxDAttr = 0;
    MassSpring::WallDataAttr *wallAttr = 0;
    MassSpring::CellModelAttr *cellDAttr = 0;
    PolarizerDiffusion *polarizerProcess = 0;

    double convergeThresh = 1.e-6;

  };


  class MassSpringSolver : public Process, public Solver<Point3d>
  {
  public:
    MassSpringSolver(const Process &process): Process(process) 
    {
      setName("Model/Cell Ovule Growth/20 Mass Spring Mechanics/f Mass Spring Solver");

      addParm("Converge Threshold", "Threshold for convergence", ".1");
      addParm("Tissue Process", "Name of the process for the Cell Tissue", "Model/Cell Ovule Growth/40 Cell Tissue/a Cell Tissue Process");
      addParm("Mass Spring Process", "Name of Mass-Spring derivatives process", "Model/Cell Ovule Growth/20 Mass Spring Mechanics/e Mass Spring Parameters");
    }
    bool initialize(QWidget *parent);
    bool step();
       
  private:
    Mesh *mesh = 0;
    CellTissueProcess *tissueProcess = 0;
    MassSpring *massSpringProcess = 0;
    //MassSpring::VertexModelAttr *vtxDAttr = 0;

    double convergeThresh = 1e-6;
  };

  class MassSpringGrowth : public Process
  {
  public:
    MassSpringGrowth(const Process &process): Process(process) 
    {
      setName("Model/Cell Ovule Growth/30 Growth Process/b Mass Spring Growth");

      addParm("Growth Signal-Based Dt", "Growth time step for signal-based growth", "0.01");
      addParm("Growth Rate Strain", "Growth rate for strain based growth", "1.0");
      addParm("Strain based growth rate par to field", "Strain based growth rate parallel to field", "0.0");
      addParm("Strain based growth rate ortho to field", "Strain based growth rate orthogonal to field", "0.0");
      addParm("Strain Threshold", "Threshold for strain for growth", ".05");
      addParm("Isotropic Growth", "Do only isotropic growth for the signal-based one", "True",  QStringList() <<"False"<< "True");
      addParm("Parallel Growth to Field", "Coefficient of growth parallel to the vector field (value between 0 and 1)", "1.");
      addParm("Orthogonal Growth to Field", "Coefficient of growth orthogonal to the vector field (value between 0 and 1)", "1.");
      addParm("Tissue Process", "Name of the process for the Cell Tissue", "Model/Cell Ovule Growth/40 Cell Tissue/a Cell Tissue Process");
    }

    bool initialize(QWidget *parent);
    bool step();

  private:
    Mesh *mesh = 0;
    CellTissueProcess *tissueProcess = 0;
    //AuxinGradient::CellDataAttr *cellADAttr = 0;
    double growthDt;
    double growthRate;
    double growthRateStrParall;
    double growthRateStrOrtho;
    double strainThresh;
    bool isoGrowth;
    double paralGrowth;
    double orthoGrowth;

  };

  /////class to draw the polarizer field obtained from the polarizer diffusion
  class PolarizerRender : public Process
  {
    public:
    enum ParmNames { pOutputCC, pDrawPolarizer, pVectorSize, pNumParms };

    PolarizerRender(const Process &process) : Process(process) 
    {
      setName("Model/Cell Ovule Growth/30 Growth Process/e Polarizer Vector Render");
      setDesc("Draw Polarizer vector from diffusion gradient for the dual graph.");
      setIcon(QIcon(":/images/Default.png"));

      //addParm("Source CC", "Name of source cell complex", "");
      addParm("Output CC", "Name of output cell complex", "Draw Polarizer");
      addParm("Draw Polarizer", "Draw polarizer vector", "Yes", booleanChoice());
      addParm("Polarizer Vector Size", "Amount to scale polarizer vector", "1.0");
    }
    bool initialize(QWidget* parent);
  
    bool run() { return run(currentMesh()); }
    bool run(Mesh *mesh);

  private:
    // Parameters
    QString SourceCC;
    QString OutputCC;
    bool DrawPolarizer;
    double PolarizerSize;
    MassSpring::CellModelAttr *cellDAttr = 0; 

  };
  ///////////////////////////////////////////////////////////////////////////////////////
  /* 
  CellTissueCelDivideOvule class
  */
  class CellTissueCellDivideOvule : public Process, virtual public Cell2dDivideParms
  {
  public:
    CellTissueCellDivideOvule(const Process &process) : Process(process) 
    {
      setName("Model/Cell Ovule Growth/40 Cell Tissue/b Divide Cells");
      setDesc("Divide cells Parameters for cell division on tissue mesh.");
      setIcon(QIcon(":/images/CellDivide.png"));
      addParm("Selective propagation", "Propagate fixed concentration for growth morphogen only if at least 2 neighbours are fixed too otherwise choose only one of daughter cells", "True", QStringList()<< "True" << "False");
      addParm("Division Algorithm", "Name of the Division Algorithm", "0", QStringList() << "0" << "1");
      addParm("Tissue Process", "Name of the Tissue process", "Model/Cell Ovule Growth/40 Cell Tissue/a Cell Tissue Process");
    }

    ~CellTissueCellDivideOvule() {}

    // Initialize simulation, called from GUI thread
    bool initialize(QWidget* parent);
  
    // Process the parameters
    bool processParms();

    // Run a step of cell division
    bool step() 
    { 
      Mesh *mesh = currentMesh();
      return step(mesh, &subdiv);
    }

    // Run a step of cell division
    virtual bool step(Mesh *mesh, Subdivide *subdiv);

    // Subdivide object
    MDXSubdivide *subdivider() { return &subdiv; }
    MassSpring::CellModelAttr *cellModelAttr = 0;
    AuxinGradient::CellDataAttr *cellDataAttr = 0;
    Point3d preferredVector;

    CellTissue *tissue;
    CellTissueProcess *cellTissueProcess;

  private: 
    //CellTissue *tissue;
    MDXSubdivide subdiv;

    // Tissue parms process
    //CellTissueProcess *cellTissueProcess;

    // Parameters
    QString TissueName;
    QString TissueDualName;
    QString TissueParmsProcessName;

    Point3d divVector;
    uint DivAlg;
    double CellMaxArea;
    double *CellMaxAreaOvule;
    //Point3d preferredVector;
    bool Verbose;
    bool selectivePropagation;
   
  };
  //////////////////////////////////////////////////////////////////////////
  


  /*
   * Cell Disk Model classes
   */
  class CellDiskTissue : public CellTissueProcess
  {
  public:
    CellDiskTissue(const Process &process) : CellTissueProcess(process)
    {
      setName("Model/Cell Ovule Growth/40 Cell Tissue/a Cell Tissue Process");
    }
  };

  class CellDiskSubdivide : virtual public mdx::Subdivide
  {
  public:
    CellDiskSubdivide() {}
    void splitCellUpdate(Dimension dim, const CCStructure &cs, const CCStructure::SplitStruct& ss, 
        CCIndex otherP = CCIndex(), CCIndex otherN = CCIndex(),  double interpPos = 0.5);

    MDXSubdivide mdx;
    MassSpring::Subdivide massSpring;
    AuxinGradient::Subdivide auxinGradient;
    CellTissue *tissue = 0;
  };

  class CellDiskDivide : public CellTissueCellDivideOvule
  {
  public:
    CellDiskDivide(const Process &process) : CellTissueCellDivideOvule(process) 
    {
      setName("Model/Cell Ovule Growth/40 Cell Tissue/b Divide Cells");

      //addParm("Preferred Direction propagation", "Preferred direction to propagate cell properties for the isolated cell", "0 1 0", QStringList() << "0 1 0" << "1 0 0" << "0 0 1" << "Random" ); 
      addParm("Tissue Process", "Name of the process for the Cell Tissue", "Model/Cell Ovule Growth/40 Cell Tissue/a Cell Tissue Process");
    }

    // Initialize to grab subdivider
    bool initialize(QWidget *parent);

    // Run a step of cell division
    bool step();
    //Point3d propagationDirection = Point3d(0., 0., 0.);
    CellDiskSubdivide &subdivider() { return subdiv; }
    //MassSpring::CellModelAttr *cellModelAttr =0;
  private:
    CellDiskSubdivide subdiv;

    CellDiskTissue *tissueProcess = 0;
  };

  class CellDiskSplitEdges : public SplitEdges
  {
  public:
    CellDiskSplitEdges(const Process &process) : SplitEdges(process) 
    {
      setName("Model/Cell Ovule Growth/40 Cell Tissue/c Split Edges");
    }
  };
  class VisualizeAssignedCellProperties : public Process
  {

    public:
      VisualizeAssignedCellProperties(const Process &process) : Process(process) 
      {
        setName("Model/Cell Ovule Growth/50 Visualize Assigned Cell Properties");
        setDesc("Visualize the different assigned parameters");

        addParm("Field to visualize", "Name of the cell model attribute to visualize", "Cell Type", QStringList() << "Cell Type" << "Max Area2Division" << "Cell Parall Growth" << "Cell Ortho Growth"<< "Cell Iso Growth"
                << "Cell Pressure" << "Edge Stiffness" << "Strain Based Growth Rate" << "Growth Rate Strain Parall" << "Growth Rate Strain Ortho" << "Strain treshold" << "Special Growth Concentration");
      }
 
      Point2d calcBoundsVis(CCIndexDoubleAttr &signalAttr);
 
      // Set the selected elements growth parameters for a transverse isotopic growth
      bool run()
      {
        Mesh *mesh = currentMesh();
        if(!mesh)
          throw QString("%1::run Invalid mesh").arg(name());
    
        // Get the attribute maps
        MassSpring::CellModelAttr &cellSpringAttr = mesh->attributes().attrMap<CCIndex, MassSpring::CellModelData>("CellModelData");
        MassSpring::VertexModelAttr &vertexSpringAttr =  mesh->attributes().attrMap<CCIndex, MassSpring::VertexModelData>("VertexModelData");

        MassSpring::WallDataAttr &wallSpringAttr =  mesh->attributes().attrMap<CCIndex, MassSpring::WallData>("MassSpringWallData");

        QString visualize = parm("Field to visualize");
        if(visualize.isEmpty())
          throw QString("%1::run Visualize choice is empty").arg(name());

        QString signal = parm("Field to visualize");
        bool result = run(*mesh, cellSpringAttr, vertexSpringAttr, wallSpringAttr, visualize, signal);
        if(result) {
          mesh->updateProperties();
          mesh->setSignal(signal);
        }
        return result;
      }
      bool run(Mesh &mesh, const MassSpring::CellModelAttr &cellSpringAttr, const MassSpring::VertexModelAttr &vertexSpringAttr, const MassSpring::WallDataAttr &wallSpringAttr, const QString &visualize, const QString &signal);
  };
  class CellDiskVisDirichlet : public fem::VisDirichlet
  {
  public:
    CellDiskVisDirichlet(const Process &proc) : VisDirichlet(proc) 
    {
      setName("Model/Cell Ovule Growth/60 Visualize Dirichlet");
    }
  };


  // Main model class
  class CellDisk : public Process
  {
  public:
    CellDisk(const Process &process) : Process(process)
    {
      setName("Model/Cell Ovule Growth/00 Ovule Master Process");

      addParm("Tissue Process", "Name of Cell Tissue Process", "Model/Cell Ovule Growth/04 Cell Tissue/a Cell Tissue Process");
      addParm("Growth Signal Process", "Name of Growth Signal Solver Process", "Model/Cell Ovule Growth/10 Growth Signal/c Growth Signal Solver");
      addParm("Set Stiffness Gradient", "Name of the Set Stiffness Process", "Model/Cell Ovule Growth/20 Mass Spring Mechanics/c Set Wall Stiffness/Assign");
      addParm("Polarizer Process", "Name of Polarizer Diffusion Solver Process", "Model/Cell Ovule Growth/30 Growth Process/d Solver Diffusion for Polarizer");
      addParm("Mass Spring Process", "Name of Mass-Spring Solver Process", "Model/Cell Ovule Growth/20 Mass Spring Mechanics/f Mass Spring Solver");
      addParm("Growth Process", "Name of the process for Growth", "Model/Cell Ovule Growth/30 Growth Process/b Mass Spring Growth");
      addParm("Divide Process", "Name of the process for Cell Division", "Model/Cell Ovule Growth/40 Cell Tissue/b Divide Cells");
      addParm("Split Edges Process", "Name of the process to split edges", "Model/Cell Ovule Growth/40 Cell Tissue/c Split Edges");
      addParm("Only L1 Periclinal Growht Process", "Name of the process that controls prevention of anticlinal growth in L1", "Model/Cell Ovule Growth/20 Mass Spring Mechanics/d Only Periclinal Growth in L1");
      addParm("Max Auxin Iterations", "Maximal number of auxin solver iteration before return", "50");
      addParm("Max Polarizer Iterations", "Maximal number of polarizer solver iteration before return", "50");
      addParm("Max Mass Spring Iterations", "Maximal number of mass spring solver iteration before return", "30");
      addParm("Simulation time to save mesh", "At which growth time points save the mesh -- need to be in crescent order", "49, 46, 32" );
      addParm("Save Mesh Process Name", "Name of the Process to save the mesh", "Mesh/System/Save" );
      addParm("Name of mesh", "Name of the mesh to be saved", "OvuleCell2DGrowth");
      addParm("Self Stop Time", "Growth time to auto-stop the simulation", "120");

    }

    // Initialize the model
    bool initialize(QWidget *parent);

    // Run a step of the model
    bool step();;

    // Rewind the model (reloads the mesh)
    bool rewind(QWidget *parent);

  private:
    Mesh *mesh = 0;
    CellDiskTissue *tissueProcess = 0;
    AuxinSolver *auxinProcess = 0;
    SetWallStiffness *setWallStiffness = 0;
    PolarizerDiffusionSolver *polarizerProcess = 0;
    MassSpringSolver *massSpringProcess = 0;
    MassSpringGrowth *growthProcess = 0;
    CellDiskDivide *divideProcess = 0;
    CellDiskSplitEdges *splitEdgesProcess = 0;
    OnlyPericlinalGrowthL1 *onlyPeriL1grProcess = 0;
    MeshSave *meshSave = 0;

    QString savingTimes;
    QStringList listSavingTimes;
    QString meshName;
    

    double maxAuxIter = 0;
    double maxPolIter = 0;
    double maxMassSpringIter = 0;
    int selfStopTime = 0;
  };
}
#endif 
