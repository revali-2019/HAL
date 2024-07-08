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
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static ABMBoneMarrow.CellType.*;
import static ABMBoneMarrow.Chemokine.CXCL9;
import static HAL.Util.*;


//////////////
//GRID CLASS//
//////////////

//Records initial conditions for many different simulations
record Experiment(int run, int activeCount, int inactiveCount, int exhaustedCount, int tceStart, int tceDuration,
                  int iterations) {
}


enum CellType {
    ACTIVE(RGB256(17, 150, 150), 0, 1.0 / 24, 0, 3.0e9, 0.1 * 3.0, 0),
    INACTIVE(RGB256(255, 165, 0), 3, 1.0 / 24, 0, 0, 0, 0),
    EXHAUSTED(RGB256(200, 50, 250), 1, 1.0 / 24, 0, 0, 0, 0),
    MYELOMA(RGB256(47, 32, 66), 2, 1.0 / 72, 1.0 / 24, 0, 0, 0.0005),
    RESISTANT(RGB256(255, 0, 0), 6, 1.0 / 72, 1.0 / 24 * 0.9, 0, 0, 0), // 90% of Myeloma division rate
    LINING(RGB256(64, 106, 151), 4, 0, 0, 0, 0, 0),
    BONE(RGB256(255, 255, 250), 5, 0, 0, 0, 0, 0),
    DEFAULT;

    int rgbColor;
    int index;
    double deathProb;
    double divProb;
    double taxisCoeff;
    double diffusionCoeff;
    double mutationProb;


    CellType(int color, int index, double deathProb, double divProb, double taxisCoeff, double diffusionCoeff, double mutationProb) {
        this.rgbColor = color;
        this.index = index;
        this.deathProb = deathProb;
        this.divProb = divProb;
        this.taxisCoeff = taxisCoeff;
        this.diffusionCoeff = diffusionCoeff;
        this.mutationProb = mutationProb; // probability that a myeloma cell will lose BCMA due to a mutation NOT CLINICALLY VERIFIED
    }

    CellType() {
        this.rgbColor = RGB256(0, 0, 0);
        this.index = -999;
        this.deathProb = 0;
        this.divProb = 0;
    }

    boolean is(CellType other) {
        return this == other;
    }


}

enum Chemokine {
    CXCL9(.1, -0.3, 2.5, 1.7),
    DEFAULT;

    double productionRate;
    double decayRate;
    double diffusionCoeff;
    double maxValue;

    Chemokine(double productionRate, double decayRate, double diffusionCoeff, double maxValue) {
        this.diffusionCoeff = diffusionCoeff;
        this.decayRate = decayRate;
        this.productionRate = productionRate;
        this.maxValue = maxValue;
    }

    Chemokine() {
    }
}

class GridLattice extends AgentGrid2D<AgentLattice> {

    int timeStep = 0; // initial timestep

    final static boolean VISUALS; // turns on/off visualizations for the model
    final static int xDim, yDim, DAYS, STEPS;
    final static List<Experiment> experimentsList;

    static {
        experimentsList = Arrays.asList(
                new Experiment(1, 25, 500, 0, 10, 15, 8)
        );
        xDim = 160;
        yDim = 150;
        DAYS = 200;
        STEPS = DAYS * 24;
        VISUALS = true;
    }


    boolean TCE_ON = true; // switching on tce therapy


    public double[] cellCounts = new double[7];
    public PDEGrid2D cxcl9; // pde grid for CXCL9

    public int[] tmoveHood = MooreHood(true); //Have option of no movement
    public Rand rn = new Rand();
    FileIO output;
    GifMaker gmCellVis, gmCXCL9Vis;
    UIGrid Cell_vis = new UIGrid(xDim, yDim, 4, 2, 5);
    UIGrid CXCL9_vis = new UIGrid(xDim, yDim, 4, 2, 5);

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

        cxcl9 = new PDEGrid2D(xDim, yDim, false, false); //This assumes PERIODIC BOUNDARY CONDITION (wrapX=TRUE,wrapY=TRUE)

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
                    vis.SetPix(x, y, HeatMapRGB(cxcl9.Get(x, y)));
                } else {
                    vis.SetPix(x, y, HeatMapRGB(cxcl9.Get(x, y)));

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
        for (int x = 0; x < cxcl9.xDim; x++) {
            for (int y = 0; y < cxcl9.yDim; y++) {
                if (GetAgent(x, y) != null && (GetAgent(x, y).cell.is(MYELOMA) || GetAgent(x, y).cell.is(RESISTANT))) {
                    cxcl9.Add(x, y, CXCL9.productionRate / CXCL9.maxValue);
                }
            }
        }


        for (int i = 0; i < xDim * yDim; i++) {
            if (this.GetAgent(i) != null) {
                this.GetAgent(i).attachTCE();
            }
        }

        //TCE Diffusion
        cxcl9.DiffusionADI(CXCL9.diffusionCoeff);


        //Natural Decay of RANKL
        cxcl9.MulAll(CXCL9.decayRate);
        cxcl9.Update();


        ShuffleAgents(rn);
        for (AgentLattice c : this) {
            c.CellStep(day);
        }
    }

    public static String makeOutputPath(Experiment e) {
        // PATHS
        String projPath = "ABMBoneMarrow";
        //Save in folder format "initActiveTCells_initInactiveTCell_initExhaustedTcell"
        String setting_dir = projPath + "/output/" + "/" + "run:" + e.run() + "act:" + e.activeCount() + "_in:" + e.inactiveCount() + "_ex:" + e.exhaustedCount() + "_start:" + e.tceStart() + "_duration:" + e.tceDuration();
        // CREATE OUTPUT DIR
        new File(setting_dir).mkdirs();
        String output_dir = setting_dir;
        new File(output_dir).mkdirs();
        return output_dir;
    }

    public static UIWindow initWindow(UIGrid Cell_vis, UIGrid CXCL9_vis, int iteration, Experiment e) {
        UIWindow win = new UIWindow("Active: " + e.activeCount() + ", Inactive: " + e.inactiveCount() + ", Exhausted: " + e.exhaustedCount() + ", Start: " + e.tceStart() + ", Duration " + e.tceDuration() + ", Iteration " + (iteration + 1));
        win.AddCol(0, new UILabel("Cells"));
        win.AddCol(0, Cell_vis);
        win.AddCol(2, new UILabel("CXCl9"));
        win.AddCol(2, CXCL9_vis);
        return win;
    }

    void generateFiles(String outputDir, int iteration) {
        output = new FileIO(outputDir.concat("/").concat("CellCounts_" + iteration).concat(".csv"), "w");
        output.Write("Time,TCELLS, EXTCELLS, MYELOMA, Inactive Tcells, CXCL9, Bone, Resistant" + "\n");

        gmCellVis = new GifMaker(outputDir.concat("/").concat("CellVid_" + iteration).concat(".gif"), 10, true);
        gmCXCL9Vis = new GifMaker(outputDir.concat("/").concat("CXCL9Vid_" + iteration).concat(".gif"), 10, true);

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

        Map<Experiment, String> experiments = experimentsList.stream().collect(Collectors.toMap(e -> e, e -> makeOutputPath(e)));

        experiments.keySet()
                .stream()
                .parallel()
                .forEach(e -> {
                    IntStream.range(0, e.iterations())
                            .parallel()
                            .forEach(i -> {
                                int day = 0;
                                boolean initial_recruitment = true;
                                var outputDir = experiments.get(e);

                                // GRID
                                GridLattice g = new GridLattice();
                                g.generateFiles(outputDir, i);

                                // OUTPUT WINDOW
                                UIWindow win = initWindow(g.Cell_vis, g.CXCL9_vis, i, e);
                                String title = "Run: " + e.run() + ", Active: " + e.activeCount() + ", Inactive: " + e.inactiveCount() + ", Exhausted: " + e.exhaustedCount() + ", Start: " + e.tceStart() + ", Duration " + e.tceDuration() + ", Iteration " + (i + 1);
                                win.RunGui();
                                // TIME LOOP
                                for (int t = 0; t < STEPS; t++) {
                                    if (g.timeStep % 48 == 0) {
                                        g.gmCellVis.AddFrame(g.Cell_vis);
                                        g.gmCXCL9Vis.AddFrame(g.CXCL9_vis);
                                    }
                                    win.frame.setTitle(title + ", Day " + day);
                                    if (VISUALS) {
                                        win.TickPause(0); //slows or speeds up the video frame rate
                                    }

                                    g.TCE_ON = (day >= e.tceStart() && day < e.tceStart() + e.tceDuration());
                                    g.timeStep = t;

                                    // COUNT CURRENT POPULATION
                                    //
                                    double[] cts = g.CellCounts();

                                    if (cts[MYELOMA.index] >= 500) { //recruiting tcells (if number of myeloma cells >= 500)
                                        if (initial_recruitment) {
                                            g.recruitCells(ACTIVE, e.activeCount());
                                        }
                                        initial_recruitment = false;
                                    }

                                    if (g.timeStep == 0) {
                                        g.recruitCells(INACTIVE, e.inactiveCount());

                                        g.recruitCells(EXHAUSTED, e.exhaustedCount());
                                    }

                                    //Putting CXCL9 in cell population map so we can use population data in CSV
                                    cts[LINING.index] = g.cxcl9.GetAvg();
                                    // RUN MODEL SIMULATION STEP
                                    g.ModelStep(g.timeStep, day);
                                    // GRAPHICAL OUTPUT
                                    if (VISUALS) {
                                        g.Draw(g.Cell_vis);
                                        g.DrawCXCL9(g.CXCL9_vis);

                                    }
                                    // OTHER OUTPUT
                                    if (g.timeStep % 24.0 == 0) {
                                        day += 1;

                                        g.RecordOut(g.output, g.timeStep / 24, cts);
                                    }
                                }

                                g.cellCounts = g.CellCounts();

                                g.output.Close();
//                                win.Close();
                                g.gmCellVis.Close();
                                g.gmCXCL9Vis.Close();

                            });

                });

//        System.exit(0);
    }

    private void recruitCells(CellType cellType, int count) {
        int k = 0;
        while (k < count) {
            int xinit = rn.Int(xDim);
            int yinit = rn.Int(yDim);
            while (this.PopAt(xinit, yinit) > 0) {
                xinit = rn.Int(xDim);
                yinit = rn.Int(yDim);
            }
            // Once xinit and yinit are within (0,0) and (xDim,yDim), place agent.
            NewAgentSQ(xinit, yinit).as(cellType);
            k++;
        }
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

    boolean TCE = false;// used for whether a tcell is attached to TCE

    ////////////////
    //CELL METHODS//
    ////////////////


    public AgentLattice as(CellType type) {
        this.cell = type;
        return this;
    }


    public int seekCXCL9() {
        int neighbors = MapHood(G.tmoveHood); // includes self
        double[] CXCL9_levels = new double[9]; // stores CXCL9 levels at the 8 directions and the center

        // Extracting CXCL9 levels for all 8 directions around the center position
        for (int i = 0; i < 9; i++) {
            CXCL9_levels[i] = G.cxcl9.Get(G.tmoveHood[i]);
        }

        double rnum = G.rn.Double();
        double ProbSum = 0.0;
        ArrayList<Integer> emptyHood = new ArrayList<>();
        ArrayList<Double> cProbArray = new ArrayList<>(); // cumulative probabilities

        for (int i = 0; i < neighbors; i++) {
            if (G.GetAgent(G.tmoveHood[i]) == null) {
                emptyHood.add(G.tmoveHood[i]); // add index to list

                double P = switch (i) {
                    case 1 -> // right
                            ACTIVE.diffusionCoeff + (ACTIVE.taxisCoeff * CXCL9.maxValue) / 8 * (CXCL9_levels[1] - CXCL9_levels[2]);
                    case 2 -> // left
                            ACTIVE.diffusionCoeff - (ACTIVE.taxisCoeff * CXCL9.maxValue) / 8 * (CXCL9_levels[1] - CXCL9_levels[2]);
                    case 3 -> // up
                            ACTIVE.diffusionCoeff + (ACTIVE.taxisCoeff * CXCL9.maxValue) / 8 * (CXCL9_levels[3] - CXCL9_levels[4]);
                    case 4 -> // down
                            ACTIVE.diffusionCoeff - (ACTIVE.taxisCoeff * CXCL9.maxValue) / 8 * (CXCL9_levels[3] - CXCL9_levels[4]);
                    case 5 -> // top-right
                            ACTIVE.diffusionCoeff + (ACTIVE.taxisCoeff * CXCL9.maxValue) / 8 * (CXCL9_levels[5] - CXCL9_levels[6]);
                    case 6 -> // top-left
                            ACTIVE.diffusionCoeff - (ACTIVE.taxisCoeff * CXCL9.maxValue) / 8 * (CXCL9_levels[5] - CXCL9_levels[6]);
                    case 7 -> // bottom-right
                            ACTIVE.diffusionCoeff + (ACTIVE.taxisCoeff * CXCL9.maxValue) / 8 * (CXCL9_levels[7] - CXCL9_levels[8]);
                    case 8 -> // bottom-left
                            ACTIVE.diffusionCoeff - (ACTIVE.taxisCoeff * CXCL9.maxValue) / 8 * (CXCL9_levels[7] - CXCL9_levels[8]);
                    default -> 0;
                };

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
        if (cell.is(MYELOMA)) {
            this.Dispose();
        }
    }

    public void attachTCE() { //Assume even distribution of TCEs
        double r = G.rn.Double();

        if (G.TCE_ON && this.cell.is(INACTIVE) && r < 0.01) {
//            G.printCellPopulation();
            this.cell = ACTIVE;
            this.TCE = true;
        }
    }


    void CellStep(double day) {

        if (cell.is(EXHAUSTED)) {

            int[] movdivHood = MooreHood(true); // for division and movement
            int emptyNeighbors = MapEmptyHood(movdivHood); // mapping empty spots
            double rn_BirthDeath = G.rn.Double();
            double pdeath = EXHAUSTED.deathProb;
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

        if (cell.is(ACTIVE)) {
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
                        if (G.GetAgent(movdivHood[j]) != null && (G.GetAgent(movdivHood[j]).cell.is(MYELOMA) || G.GetAgent(movdivHood[j]).cell.is(RESISTANT)) && this.hours_since_kill >= 1) {
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


        if (cell.is(MYELOMA)) {
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
                    child.as(MYELOMA);
                    if (G.rn.Double() < MYELOMA.mutationProb) {
                        child.as(RESISTANT);
//                    } else if (G.rn.Double() > mutProb) {
                    }
                }
            }
        }

        if (cell.is(RESISTANT)) {
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

        if (cell.is(INACTIVE)) {
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


