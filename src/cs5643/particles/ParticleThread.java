package cs5643.particles;

import com.sun.org.apache.xml.internal.dtm.ref.DTMDefaultBaseIterators;
import jdk.internal.org.objectweb.asm.tree.TryCatchBlockNode;

import javax.media.opengl.GL2;
import java.security.cert.TrustAnchor;
import java.util.concurrent.CountDownLatch;

public class ParticleThread implements Runnable {


    public static CountDownLatch latch;
    public static Task task;
	public static ParticleSystem particleSystem;
	public static GL2 gl;
    private Particle p;

    ParticleThread(Particle p){
	    this.p = p;
    }

    public void run(){
	    try {
		    switch (task) {
			    case UpdateVandXStar:
				    particleSystem.updateVandXstar(p);
				    break;
			    case UpdateCells:
				    particleSystem.updateCells(p);
				    break;
			    case CalculateLambda:
				    particleSystem.calculateLambda(p);
				    break;
			    case CalculateDp:
				    particleSystem.calculateDp(p);
				    break;
			    case XSPHViscosity:
				    particleSystem.XSPHViscosity(p);
				    break;
			    case VorticityConfinement:
				    particleSystem.VorticityConfinement(p);
				    break;
			    case Display:
				    p.display(gl);
				    break;
		    }
	    } catch (InterruptedException e){
		    e.printStackTrace();
	    }
	    latch.countDown();
    }

	public static void init(ParticleSystem ps) {
		particleSystem = ps;
		task = Task.None;
	}



    public enum Task {
	    None,
        UpdateVandXStar,
        UpdateCells,
        CalculateLambda,
        CalculateDp,
        XSPHViscosity,
        VorticityConfinement,
	    Display
    }
}
