package cs5643.particles;

import org.omg.PortableInterceptor.SYSTEM_EXCEPTION;

import javax.vecmath.*;
import javax.media.opengl.*;
import javax.media.opengl.glu.*;

/**
 * Simple particle implementation, with miscellaneous adornments.
 *
 * @author Doug James, January 2007
 * @author Eston Schweickart, February 2014
 */
public class Particle
{
    /** Radius of particle's sphere graphic. */
    public static final double PARTICLE_RADIUS = 0.015;

    /** Display list index. */
    private static int PARTICLE_DISPLAY_LIST = -1;

    /** Highlighted appearance if true, otherwise white. */
    private boolean highlight = false;

    /** Default mass. */
    double   m = Constants.PARTICLE_MASS;

    /** Deformed Position. */
    Point3d  x = new Point3d();

    /** Undeformed/material Position. */
    Point3d  x0 = new Point3d();

	/** Density at this particle */
	double rho = 0;

	/** Temporary position for simulation loop */
	Point3d xStar = new Point3d();

	/** Change in position for simulation loop */
	Vector3d dP = new Vector3d();

    /** Velocity. */
    Vector3d v = new Vector3d();

	/** Proposed velocity */
	Vector3d vStar = new Vector3d();

	/** Vorticity */
	Vector3d omega = new Vector3d();

    /** Force accumulator. */
    Vector3d f = new Vector3d();

	/** Color of particle */
	public float[] c;

	/** Cell this particle belongs to */
	public int[] cell;

	/** Lambda of each particle */
	public double lambda;

    /** ParticleThread for threading. Won't happen in future iteration :P */
    public ParticleThread thread;

    /**
     * Constructs particle with the specified material/undeformed
     * coordinate, x0.
     */
    Particle(Point3d x0)
    {
        this.x0.set(x0);
        x.set(x0);
	    c = new float[4];
	    c[0] = 0; c[1] = 1; c[2] = 1; c[3] = 1;
	    cell = new int[3];
        thread = new ParticleThread(this);
    }

    /** Draws spherical particle using a display list. */
    public void display(GL2 gl)
    {
        if(PARTICLE_DISPLAY_LIST < 0) {// MAKE DISPLAY LIST:
            int displayListIndex = gl.glGenLists(1);
            GLU glu = GLU.createGLU();
            GLUquadric quadric = glu.gluNewQuadric();
            gl.glNewList(displayListIndex, GL2.GL_COMPILE);
            glu.gluSphere(quadric, PARTICLE_RADIUS, 16, 8);
            gl.glEndList();
            glu.gluDeleteQuadric(quadric);

            // For older versions of JOGL
            //glu.destroy();

            System.out.println("MADE DISPLAY LIST "+displayListIndex+" : "+gl.glIsList(displayListIndex));
            PARTICLE_DISPLAY_LIST = displayListIndex;
        }

        /// COLOR: DEFAULT CYAN; GREEN IF HIGHLIGHTED
        if(highlight) {
            c[2] = 0;
        } else {
	        c[2] = 1;
        }

        // Hack to make things more colorful/interesting
        c[1] = (float)x.y;
        gl.glMaterialfv(GL2.GL_FRONT, GL2.GL_DIFFUSE, c, 0); // Color used by shader

        /// DRAW ORIGIN-CIRCLE TRANSLATED TO "p":
        gl.glPushMatrix();
        gl.glTranslated(x.x, x.y, x.z);
        gl.glCallList(PARTICLE_DISPLAY_LIST); // Draw the particle
        gl.glPopMatrix();
    }

    /** Specifies whether particle should be drawn highlighted. */
    public void setHighlight(boolean highlight) {
        this.highlight = highlight;
    }
    /** True if particle should be drawn highlighted. */
    public boolean getHighlight() {
        return highlight;
    }
}
