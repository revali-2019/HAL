package HIP_IMO;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.Gui.GridWindow;
import HAL.Rand;
import HAL.Util;

import java.util.Random;

class Cell extends AgentSQ2Dunstackable<ExampleGrid> {


    public void stepCell(double divProb, double dieProb) {
        //To access things in GRid class, use G.
        //Event: Die
        if (G.rng.Double() < dieProb) {
            Dispose();
            return;
        }

//        if (G.rng.Double() < mutateProb) {
//
//        }

        //Event: Divide
        if (G.rng.Double() < divProb) {
            int options = MapEmptyHood(G.divHood);
            if (options > 0) {
                G.NewAgentSQ(G.divHood[G.rng.Int(options)]);
            }
        }


    }


}

public class ExampleGrid extends AgentGrid2D<Cell> {
    Rand rng = new Rand();
    Random rand = new Random();

    int[] divHood = Util.VonNeumannHood(false);

    public void step(double divProb, double dieProb) {
        for (Cell c : this) {
            c.stepCell(divProb, dieProb);
        }
    }

    public void draw(GridWindow win) {
        for (int i = 0; i < length; i++) {
            int color = Util.BLACK;
            if (GetAgent(i) != null) {
                color = Util.WHITE;
            }
            win.SetPix(i, color);
        }
    }

    public static void main(String[] args) {
        ExampleGrid g = new ExampleGrid(100, 100);
        int t = 0, tMax = 1000;
        double divProb = 0.5, dieProb = 0.1;
        GridWindow win = new GridWindow(100, 100, 4);

        g.NewAgentSQ(g.xDim / 2, g.yDim / 2);
        while (t < tMax) {
            win.TickPause(50);
            g.draw(win);
            g.step(divProb, dieProb);
            t++;
        }
    }

    public ExampleGrid(int x, int y) {
        super(x, y, Cell.class);
    }
}
