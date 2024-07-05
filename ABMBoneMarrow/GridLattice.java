package ABMBoneMarrow;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GifMaker;
import HAL.Gui.UIGrid;
import HAL.Gui.UILabel;
import HAL.Gui.UIWindow;
import HAL.Rand;
import HAL.Tools.FileIO;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

import static ABMBoneMarrow.CellType.*;
import static ABMBoneMarrow.GridLattice.mutProb;
import static HAL.Util.*;

//////////////
//GRID CLASS//
//////////////

//Records initial conditions for many different simulations
record Experiments(int activeCount, int inactiveCount, int exhaustedCount, int start, int duration, int iterations) {
}


enum CellType {
    ACTIVE(RGB256(17, 150, 150), 0, 1.0 / 24, 0),
    INACTIVE(RGB256(255, 165, 0), 3, 1.0 / 24, 0),
    EXHAUSTED(RGB256(200, 50, 250), 1, 1.0 / 24, 0),
    MYELOMA(RGB256(47, 32, 66), 2, 1.0 / 72, 1.0 / 24),
    RESISTANT(RGB256(255, 0, 0), 6, 1.0 / 72, 1.0 / 24 * 0.9), // 90% of Myeloma division rate
    LINING(RGB256(64, 106, 151), 4, 0, 0),
    BONE(RGB256(255, 255, 250), 5, 0, 0),
    DEFAULT;

    int rgbColor;
    int index;

    double deathProb;

    double divProb;


    CellType(int color, int index, double deathProb, double divProb) {
        this.rgbColor = color;
        this.index = index;
        this.deathProb = deathProb;
        this.divProb = divProb;
    }

    CellType() {
        this.rgbColor = RGB256(0, 0, 0);
        this.index = -999;
        this.deathProb = 0;
        this.divProb = 0;
    }
}

class GridLattice extends AgentGrid2D<AgentLattice> {

    static final int xDim, yDim, DAYS, STEPS;
    int timeStep = 0; // initial timestep

    public final static boolean VISUALS; // turns on/off visualizations for the model

    static {
        experimentsList = Arrays.asList(
                new Experiments(1, 2, 3, 4, 5, 5),
                new Experiments(1, 2, 3, 4, 5, 5),
                new Experiments(1, 2, 3, 4, 5, 5),
                new Experiments(1, 2, 3, 4, 5, 5),
                new Experiments(1, 2, 3, 4, 5, 5)
        );
        xDim = 160;
        yDim = 150;
        DAYS = 450;
        STEPS = DAYS * 24;
        VISUALS = true;
    }
    //TODO THE ACTIVATION PROBLEM HAS SOMETHING TO DO WITH TIMESTEP BECOMING 1

    public final static int ITERATIONS = 10; // number of times the model will run
    public static final int InitactiveTcells = 0; // initial number of active tcells
    public static final int InitinactiveTcells = 500; //initial number of inactive tcells
    public static final int InitexhaustedTcells = 0; // initial number of exhausted tcells
    final static int TCE_start = 5; // start day for treatment

    final static int TCE_Duration = 50; // duration (in days) of treatment


    //TODO Fix resistant growth rate to be smaller,
    final static double mutProb = 0.0005; // probability that a myeloma cell will lose BCMA due to a mutation NOT CLINICALLY VERIFIED
    public double T_CELL_DEATH_RATE = 1.0 / 24; // tcell death rate - used for exhausted tcells


    double CXCL9_productionRate = .1; // production of CXCL9
    double CXCL9_decayRate = -0.3; // decay of CXCL9
    double CXCL9_DiffCoef = 2.5; // diffusion coefficient for CXCL9
    double maxCXCL9 = 1.7; // used for numerical stability
    double Tcell_TaxisCoeff = 3.0e9; // tcell movement towards CXCL9

    double Tcell_DiffCoef = 0.1 * 3.0; // diffusion coefficient for cd8 t-cells

    boolean TCE_ON = true; // switching on tce therapy

    boolean preRes = true;

    double resFrac = 0;


    public double[] cellCounts = new double[7];
    public PDEGrid2D CXCL9; // pde grid for CXCL9

    public PDEGrid2D TCE; // TCE distribution pde grid

    public int[] tmoveHood = MooreHood(true); //Have option of no movement


    public Rand rn = new Rand();

    public static List<Experiments> experimentsList;


    FileIO output;
    ////////////////////
    //GRID CONSTRUCTOR//

    ////////////////////
    //FIXME play around with bounds T Cell kill rate

    public double boundedGaussian(double mean, double dev, double min, double max) {
        double gauss = rn.Gaussian(0, 1);
        double val = dev * gauss + mean;
        while (val > max || val < min) {
            gauss = rn.Gaussian(0, 1);
            val = dev * gauss + mean;
        }
        return val;
    }

    boolean isWithinDistanceFromBone(int x, int y, int startX, int startY, int sideLength, int distance) {
        int boneEndX = startX + sideLength;
        int boneEndY = startY + sideLength;

        int minX = Math.max(0, startX - distance);
        int maxX = Math.min(xDim - 1, boneEndX + distance);
        int minY = Math.max(0, startY - distance);
        int maxY = Math.min(yDim - 1, boneEndY + distance);

        return x >= minX && x <= maxX && y >= minY && y <= maxY;
    }

    public GridLattice() {

        super(GridLattice.xDim, GridLattice.yDim, AgentLattice.class, true, true);

        CXCL9 = new PDEGrid2D(xDim, yDim, false, false); //This assumes PERIODIC BOUNDARY CONDITION (wrapX=TRUE,wrapY=TRUE)

        double totalArea = xDim * yDim;
        double blackBoxArea = 0.129 * totalArea;
        int sideLength = (int) Math.sqrt(blackBoxArea);

        // Calculate the center of the grid
        int centerX = xDim / 2;
        int centerY = yDim / 2;

        // Calculate the top-left corner of the black box
        int startX = centerX - sideLength / 2;
        int startY = centerY - sideLength / 2;

        instantiateBoneCells(startX, startY, sideLength);
        instantiateMyelomaCells(startX, startY, sideLength);

        this.timeStep = 0;
    }

    private void instantiateMyelomaCells(int startX, int startY, int sideLength) {
        for (int i = 0; i < 500; i++) {
            // Recruit one myeloma cell within 10 grid spaces from the bone
            int myelomaX, myelomaY;
            do {
                myelomaX = rn.Int(xDim);
                myelomaY = rn.Int(yDim);
            } while (!isWithinDistanceFromBone(myelomaX, myelomaY, startX, startY, sideLength, 5) || PopAt(myelomaX, myelomaY) > 0);

            AgentLattice c = NewAgentSQ(myelomaX, myelomaY).as(MYELOMA);

        }
    }

    private void instantiateBoneCells(int startX, int startY, int sideLength) {
        for (int i = startX; i < startX + sideLength; i++) {
            for (int j = startY; j < startY + sideLength; j++) {
                if (i == startX || i == startX + sideLength - 1 || j == startY || j == startY + sideLength - 1) {
                    // Set the cell type to BLACK_BOX and color to black
                    NewAgentSQ(i, j).as(CellType.LINING);
                } else {
                    NewAgentSQ(i, j).as(BONE);
                }
            }
        }
    }
    ////////////////
    //GRID METHODS//
    ////////////////

    public void Draw(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                AgentLattice drawMe = GetAgent(x, y);
                if (drawMe != null && drawMe.cell != null) {
                    vis.SetPix(x, y, drawMe.cell.rgbColor);
                } else {
                    vis.SetPix(x, y, RGB256(240, 220, 220));
                }
            }
        }
    }

    public void DrawCXCL9(UIGrid vis) {
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                AgentLattice drawMe = GetAgent(x, y);
                if (drawMe != null) {
                    vis.SetPix(x, y, HeatMapRGB(CXCL9.Get(x, y)));
                } else {
                    vis.SetPix(x, y, HeatMapRGB(CXCL9.Get(x, y)));

                }
            }
        }
    }


    public double[] CellCounts() {
        double[] cts = new double[7];
        for (AgentLattice c : this) {
            if (c.cell != null) {
                cts[c.cell.index]++;
            }
        }

        return cts;
    }

    public void RecordOut(FileIO writeHere, int time, double[] cts) {
        writeHere.Write(time + ",");
        writeHere.WriteDelimit(cts, ",");
        writeHere.Write("," + "\n");
    }


    public void ModelStep(double timeStep, int day) {
        for (int x = 0; x < CXCL9.xDim; x++) {
            for (int y = 0; y < CXCL9.yDim; y++) {
                if (GetAgent(x, y) != null && (GetAgent(x, y).cell == MYELOMA || GetAgent(x, y).cell == RESISTANT)) {
                    CXCL9.Add(x, y, CXCL9_productionRate / maxCXCL9);

                }

            }
        }


        for (int i = 0; i < xDim * yDim; i++) {
            if (this.GetAgent(i) != null) {
                this.GetAgent(i).TCE_Attach();
            }
        }

        //TCE Diffusion
        CXCL9.DiffusionADI(CXCL9_DiffCoef);


        //Natural Decay of RANKL
        CXCL9.MulAll(CXCL9_decayRate);
        CXCL9.Update();


        ShuffleAgents(rn);
        for (AgentLattice c : this) {
            c.CellStep(day);
        }
    }

    public static String makeOutputPath() {
        // PATHS
        String projPath = "ABMBoneMarrow";
        //Save in folder format "initActiveTCells_initInactiveTCell_initExhaustedTcell"
        String setting_dir = projPath + "/output/" + "/" + "act:" + InitactiveTcells + "_in:" + InitinactiveTcells + "_ex:" + InitexhaustedTcells + "_start:" + (int) TCE_start + "_duration:" + (int) TCE_Duration;
        // CREATE OUTPUT DIR
        new File(setting_dir).mkdirs();
        String output_dir = setting_dir;
        new File(output_dir).mkdirs();

        return output_dir;
    }

    public static UIWindow initWindow(UIGrid Cell_vis, UIGrid CXCL9_vis, int j) {
        UIWindow win = new UIWindow("Active: " + InitactiveTcells + ", Inactive: " + InitinactiveTcells + ", Exhausted: " + InitexhaustedTcells + ", Start: " + TCE_start + ", Duration " + TCE_Duration + ", Iteration " + (j + 1));
        win.AddCol(0, new UILabel("Cells"));
        win.AddCol(0, Cell_vis);
        win.AddCol(2, new UILabel("CXCl9"));
        win.AddCol(2, CXCL9_vis);
        return win;
    }

    void generateFile(String fileWithFullPath) {
        output = new FileIO(fileWithFullPath, "w");
        output.Write("Time,TCELLS, EXTCELLS, MYELOMA, Inactive Tcells, CXCL9, Bone, Resistant" + "\n");
    }


    void printCellPopulation() {
        Map<CellType, Integer> cells = new HashMap<>();
        for (CellType type : CellType.values()) {
            cells.put(type, 0);
        }
        for (AgentLattice a : this) {

            cells.put(a.cell, cells.get(a.cell) + 1);

        }
        System.out.println(cells);
        System.out.println(timeStep);
    }

    /**
     * Visualization of tumor microenvironment for a certain number of iterations.
     * 1. Active T Cells
     * 2. Multiple Myeloma Cells
     * 3. Inactive T Cells
     * 4. Exhausted T Cells
     * 5. CXCL9 Chemokine
     *
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

        String output_dir = makeOutputPath();
        String path_to_output_file = "";


        IntStream.range(0, ITERATIONS)
                .parallel()
                .forEach(j -> {

                            int day = 0;
                            boolean initial_recruitment = true;

                            // GRID
                            GridLattice g = new GridLattice();
                            g.generateFile(output_dir.concat("/").concat("CellCounts_" + j).concat(".csv"));

                            GifMaker gm_Cell_vis = new GifMaker(output_dir.concat("/").concat("CellVid_" + j).concat(".gif"), 10, true);
                            GifMaker gm_CXCL9_vis = new GifMaker(output_dir.concat("/").concat("CXCL9Vid_" + j).concat(".gif"), 10, true);

                            // OUTPUT WINDOW
                            UIGrid Cell_vis = new UIGrid(xDim, yDim, 4, 2, 5);
                            UIGrid CXCL9_vis = new UIGrid(xDim, yDim, 4, 2, 5);
                            UIWindow win = initWindow(Cell_vis, CXCL9_vis, j);
                            String title = "Active: " + InitactiveTcells + ", Inactive: " + InitinactiveTcells + ", Exhausted: " + InitexhaustedTcells + ", Start: " + TCE_start + ", Duration " + TCE_Duration + ", Iteration " + (j + 1);
                            win.RunGui();


                            // TIME LOOP
                            for (int i = 0; i < STEPS; i++) {
                                if (i % 48 == 0) {
                                    gm_Cell_vis.AddFrame(Cell_vis);
                                    gm_CXCL9_vis.AddFrame(CXCL9_vis);
                                }
                                win.frame.setTitle(title + ", Day " + day);
//
                                if (day >= TCE_start && day < TCE_start + TCE_Duration) { //Turn on TCE treatment if day is in treatment block
                                    g.TCE_ON = true;

                                } else {
                                    g.TCE_ON = false;
                                }

                                g.timeStep = i;

                                if (VISUALS) {
                                    win.TickPause(10); //slows or speeds up the video frame rate
                                }
                                // COUNT CURRENT POPULATION
                                //
                                double[] cts = g.CellCounts();
                                //
                                if (cts[MYELOMA.index] >= 500) { //recruiting tcells (if number of myeloma cells >= 500)
                                    if (initial_recruitment) {
                                        int InitTcells = g.InitactiveTcells;
                                        int k = 0;
                                        while (k < InitTcells) {
                                            int xinit = g.rn.Int(xDim);
                                            int yinit = g.rn.Int(yDim);
                                            while (g.PopAt(xinit, yinit) > 0) {
                                                xinit = g.rn.Int(xDim);
                                                yinit = g.rn.Int(yDim);
                                            }
                                            // Once xinit and yinit are within (0,0) and (xDim,yDim), place agent.
                                            g.NewAgentSQ(xinit, yinit).as(ACTIVE);
                                            k++;
                                        }
                                    }
                                    initial_recruitment = false;
                                }

                                //Initializing Inactive TCElls
                                if (g.timeStep == 0) {
                                    int k = 0;
                                    while (k < g.InitinactiveTcells) {
                                        int xinit = g.rn.Int(xDim);
                                        int yinit = g.rn.Int(yDim);
                                        while (g.PopAt(xinit, yinit) > 0) {
                                            xinit = g.rn.Int(xDim);
                                            yinit = g.rn.Int(yDim);
                                        }
                                        // Once xinit and yinit are within (0,0) and (xDim,yDim), place agent.
                                        g.NewAgentSQ(xinit, yinit).as(INACTIVE);
                                        k++;
                                    }

                                    //Initializing Exhausted TCells
                                    int m = 0;
                                    while (m < g.InitexhaustedTcells) {
                                        int xinit = g.rn.Int(xDim);
                                        int yinit = g.rn.Int(yDim);
                                        while (g.PopAt(xinit, yinit) > 0) {
                                            xinit = g.rn.Int(xDim);
                                            yinit = g.rn.Int(yDim);
                                        }
                                        // Once xinit and yinit are within (0,0) and (xDim,yDim), place agent.
                                        g.NewAgentSQ(xinit, yinit).as(EXHAUSTED);
                                        m++;
                                    }

                                }

                                //Putting CXCL9 in cell population map so we can use population data in CSV
                                cts[LINING.index] = g.CXCL9.GetAvg();
                                // RUN MODEL SIMULATION STEP
                                g.ModelStep(g.timeStep, day);
                                // GRAPHICAL OUTPUT
                                if (VISUALS) {
                                    g.Draw(Cell_vis);
                                    g.DrawCXCL9(CXCL9_vis);

                                }
                                // OTHER OUTPUT
                                if (i % 24.0 == 0) {
                                    day += 1;

                                    g.RecordOut(g.output, i / 24, cts);
                                }
                            }


//                            System.out.println("Iteration " + j);
                            g.cellCounts = g.CellCounts();

                            g.output.Close();

                        }
                );

    }
}

//////////////
//CELL CLASS//
//////////////

class AgentLattice extends AgentSQ2Dunstackable<GridLattice> {
    public CellType cell; // type of cell
    double pd_1 = 0; // pd_1 for t cells
    double pd_l1 = 0; // pd_l1 for t cells
    double tcellAge = 0; // used for tracking the age of a tcell
    double hours_since_kill = 0; // hours since a tcell last killed a myeloma cell
    double hours_since_division = 0; //days since the Tcell has divided
    boolean BCMA = false; // used for whether or not a myeloma cell is expressing bcma

    public boolean TCE = false;// used for whether a tcell is attached to TCE

    ////////////////
    //CELL METHODS//
    ////////////////

    public boolean BCMA_NOT_SET = false;

    public AgentLattice as(CellType type) {
        this.cell = type;
        return this;
    }


    public int seekCXCL9() {
        int neighbors = MapHood(G.tmoveHood); // includes self
        double[] CXCL9_levels = new double[9]; // stores CXCL9 levels at the 8 directions and the center

        // Extracting CXCL9 levels for all 8 directions around the center position
        for (int i = 0; i < 9; i++) {
            CXCL9_levels[i] = G.CXCL9.Get(G.tmoveHood[i]);
        }

        double rnum = G.rn.Double();
        double ProbSum = 0.0;
        ArrayList<Integer> emptyHood = new ArrayList<>();
        ArrayList<Double> cProbArray = new ArrayList<>(); // cumulative probabilities

        for (int i = 0; i < neighbors; i++) {
            if (G.GetAgent(G.tmoveHood[i]) == null) {
                emptyHood.add(G.tmoveHood[i]); // add index to list

                double P = 0;
                switch (i) {
                    case 1: // right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[1] - CXCL9_levels[2]);
                        break;
                    case 2: // left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[1] - CXCL9_levels[2]);
                        break;
                    case 3: // up
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[3] - CXCL9_levels[4]);
                        break;
                    case 4: // down
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[3] - CXCL9_levels[4]);
                        break;
                    case 5: // top-right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[5] - CXCL9_levels[6]);
                        break;
                    case 6: // top-left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[5] - CXCL9_levels[6]);
                        break;
                    case 7: // bottom-right
                        P = G.Tcell_DiffCoef + (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[7] - CXCL9_levels[8]);
                        break;
                    case 8: // bottom-left
                        P = G.Tcell_DiffCoef - (G.Tcell_TaxisCoeff * G.maxCXCL9) / 8 * (CXCL9_levels[7] - CXCL9_levels[8]);
                        break;
                }

                P = Math.max(P, 0); // ensure non-negative probability

                if (cProbArray.size() > 0) {
                    cProbArray.add(P + cProbArray.get(cProbArray.size() - 1));
                } else {
                    cProbArray.add(P);
                }
                ProbSum += P;
            }
        }

        int moveToIndex = G.tmoveHood[G.rn.Int(neighbors)]; // Default to the center index if no movement

        for (int i = 0; i < cProbArray.size(); i++) {
            if (rnum <= cProbArray.get(i) / ProbSum) {
                moveToIndex = emptyHood.get(i);
                break;
            }
        }

        return moveToIndex;
    }


    void Tcell_Kill() {
        if (cell == MYELOMA) {
            this.Dispose();
        }
    }

    public void TCE_Attach() { //Assume even distribution of TCEs
        double r = G.rn.Double();

        if (G.TCE_ON && this.cell == INACTIVE && r < 0.01) {
//            G.printCellPopulation();
            this.cell = ACTIVE;
            this.TCE = true;
        }
    }


    void CellStep(double day) {

        if (cell == EXHAUSTED) {

            int[] movdivHood = MooreHood(true); // for division and movement
            int emptyNeighbors = MapEmptyHood(movdivHood); // mapping empty spots
            double rn_BirthDeath = G.rn.Double();
            double pdeath = G.T_CELL_DEATH_RATE;
            if (emptyNeighbors > 0) {
                int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                int chosenCell = movdivHood[chosenIndex]; // Get the chosen empty cell
                if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                    MoveSQ(chosenCell);
                }
            }
            if (rn_BirthDeath < ProbScale(pdeath, 1.0 / 24)) {
                this.Dispose();
            }
        }

        if (cell == ACTIVE) {
            boolean encounteredMyeloma = false; // Track if a myeloma cell is encountered
            boolean agedDeath = false; // Track if the T cell is aged and should die

            // Initialize pd_1 value if the T cell is newly created
            if (this.tcellAge == 0) {
                this.pd_1 = G.boundedGaussian(10, 5, 5, 20);
            }


            // Increment T cell age and hours since last kill every hour
            this.tcellAge += 1;
            this.hours_since_kill += 1;
            this.hours_since_division += 1;

            for (int run = 0; run < 3; run++) {
                int[] movdivHood = MooreHood(true); // For division and movement
                int options = MapOccupiedHood(movdivHood); // Mapping occupied spots
                int emptyNeighbors = MapEmptyHood(movdivHood); // Mapping empty spots

                // Age-related death logic
                if (this.tcellAge == 720 && emptyNeighbors > 0) { // Assuming 720 hours is the lifespan (30 days)
                    for (int i = 0; i < emptyNeighbors; i++) {
                        int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                        int chosenCell = movdivHood[chosenIndex]; // Get the chosen empty cell
                        if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                            AgentLattice child = G.NewAgentSQ(chosenCell).as(this.cell);


                            child.tcellAge = 0;
                            child.pd_1 = Math.floor(this.pd_1 / 4);
                            child.pd_l1 = this.pd_l1;
                            child.hours_since_division = 0;
                            break;
                        }
                    }
                    this.Dispose();
                    break;
                }

                // Transition to exhausted T cell if pd_l1 exceeds pd_1
                if (this.pd_l1 > this.pd_1) {
                    this.cell = EXHAUSTED;

                    break;
                } else {
                    // Killing logic
                    for (int j = 0; j < options; j++) {
                        if (G.GetAgent(movdivHood[j]) != null && (G.GetAgent(movdivHood[j]).cell == MYELOMA || G.GetAgent(movdivHood[j]).cell == RESISTANT) && this.hours_since_kill >= 1) {
                            if (this.TCE && !G.GetAgent(movdivHood[j]).BCMA) {
                                ;
                            } else { // If the cells don't have TCE or have BCMA
                                G.GetAgent(movdivHood[j]).Tcell_Kill();

                                this.pd_l1 += 0.5;
                                this.hours_since_kill = 0; // Reset kill timer
                            }


                            // Division logic
                            if (this.hours_since_division >= 24 && emptyNeighbors > 0) {
                                for (int i = 0; i < emptyNeighbors; i++) {
                                    int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                                    int chosenCell = movdivHood[chosenIndex]; // Get the chosen empty cell
                                    if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                                        AgentLattice child = G.NewAgentSQ(chosenCell);
                                        child.cell = this.cell;

                                        child.tcellAge = 0;
                                        child.pd_1 = Math.floor(this.pd_1 / 4);
                                        child.pd_l1 = this.pd_l1;
                                        child.hours_since_division = 0;
                                        this.hours_since_division = 0;
                                        break;
                                    }
                                }
                            }
                            encounteredMyeloma = true; // Set the flag to true if a myeloma cell is encountered
                            break; // Exit the inner loop if a myeloma cell is encountered
                        }
                    }
                }

                // If a myeloma cell is encountered, break out of the outer loop
                if (encounteredMyeloma) {
                    break;
                }

                // Movement logic
                int moveToIndex = seekCXCL9(); // Explicitly define since seekTGFB has random values
                if (G.GetAgent(moveToIndex) == null) {
                    MoveSQ(moveToIndex);
                }
            }
        }


        if (cell == MYELOMA) {
            this.BCMA = true;

            //STEP 1: death
            if (G.rn.Double() < ProbScale(MYELOMA.deathProb, 1)) {
                Dispose();
                return;
            }


            //STEP 2: division
            if (G.rn.Double() < ProbScale(MYELOMA.divProb, 1)) {
                int[] divHood = MooreHood(true);
                int emptyNeighbors = MapEmptyHood(divHood);
                //Create new Agent
                if (emptyNeighbors > 0) {
                    // TODO Filter usage
                    AgentLattice child = G.NewAgentSQ(divHood[G.rn.Int(emptyNeighbors)]);
                    if (G.rn.Double() < mutProb) {
                        child.cell = RESISTANT;
                    } else if (G.rn.Double() > mutProb) {
                        child.cell = MYELOMA;
                    }


                }
            }
        }

        if (cell == RESISTANT) {
            this.BCMA = false;

            //STEP 1: death
            if (G.rn.Double() < ProbScale(RESISTANT.deathProb, 1)) {
                Dispose();
                return;
            }


            //STEP 2: division
            if (G.rn.Double() < ProbScale(RESISTANT.divProb, 1)) {
                int[] divHood = MooreHood(true);
                int emptyNeighbors = MapEmptyHood(divHood);
                //Create new Agent
                if (emptyNeighbors > 0) {
                    AgentLattice child = G.NewAgentSQ(divHood[G.rn.Int(emptyNeighbors)]);
                    child.cell = this.cell;

                }
            }

        }

        if (cell == INACTIVE) {
            boolean agedDeath = false; //Tracking for aged Tcells
            if (G.timeStep % 24.0 == 0) {
                this.tcellAge += 1;
            }


            for (int run = 0; run < 3; run++) {
                double rn_BirthDeath = G.rn.Double();
                //double pdiv = G.T_CELL_DIV_RATE;
                int[] movdivHood = MooreHood(true); // for division and movement
                int options = MapOccupiedHood(movdivHood); // mapping occupied spots
                int emptyNeighbors = MapEmptyHood(movdivHood); // mapping empty spots
                if (this.tcellAge == 720) {
                    if (emptyNeighbors > 0) {
                        for (int i = 0; i < emptyNeighbors; i++) {
                            int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                            int chosenCell = movdivHood[chosenIndex]; // Get the chosen empty cell
                            if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                                AgentLattice child = G.NewAgentSQ(chosenCell);
                                child.cell = this.cell;
                                child.tcellAge = 0;
                                child.pd_1 = 0;
                                child.pd_l1 = 0;
                                child.hours_since_division = 0;

                                break;
                            }
                        }
                        agedDeath = true;
                        this.Dispose();
                    }
                }
                if (agedDeath) {
                    break;
                }
                if (emptyNeighbors > 0) {
                    int chosenIndex = G.rn.Int(emptyNeighbors); // Randomly choose an empty cell index
                    int chosenCell = movdivHood[chosenIndex]; // Get the chosen empty cell
                    if (G.GetAgent(chosenCell) == null) { // Check if the chosen cell is still empty
                        MoveSQ(chosenCell);
                    }
                }
            }
        }
    }

}


