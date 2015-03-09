package cs5643.particles;

import java.util.concurrent.CountDownLatch;

public class ParticleThread implements Runnable {


    public static CountDownLatch latch;
    private static Task task;
    private Particle p;

    ParticleThread(Particle p){
        this.p = p;
    }

    public void run(){
        switch (task) {
            case UpdateVandXStar:
                break;
            case UpdateCells:
                break;
            case CalculateLambda:
                break;
            case CalculateDp:
                 break;
            case XSPHViscosity:
                break;
            case VorticityConfinement:
                break;
        }
    }

    public enum Task {
        UpdateVandXStar,
        UpdateCells,
        CalculateLambda,
        CalculateDp,
        XSPHViscosity,
        VorticityConfinement
    }
}
