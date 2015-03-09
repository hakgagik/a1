package cs5643.particles;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import javax.vecmath.*;
import javax.media.opengl.*;
import com.jogamp.opengl.util.glsl.*;


/**
 * Maintains dynamic lists of Particle and Force objects, and provides
 * access to their state for numerical integration of dynamics.
 *
 * @author Doug James, January 2007
 * @author Eston Schweickart, February 2014
 */
public class ParticleSystem //implements Serializable
{
    /** Current simulation time. */
    public static double time = 0;


	/** Number of iterations of the main iteration loop */
	int maxIter = 4;

	/** c used for XPHS viscosity */
	public final double c = 0.05;

	/** Rest density */
	public final double rho_0 = 6378.0;

	/** CFM Parameter */
	public final double epsilon = 600;

	/** Vorticity epsilon */
	public final double eVorticity = 2;

	/** Artificial pressure terms */
	public final double k = -0.0001;
	public final double q;

	/** Dimensions of the simulation cell. */
	public static final double dims = 4;
    public static final double bounds1 = dims * 0.4;
    public static final double bounds2 = dims * 0.8;

    /** Number of processors available */
    public static final int numProcessors = Runtime.getRuntime().availableProcessors();

	/** Strength of the whirlpool */
	public static final double whirl = 5;

    /** List of Particle objects. */
    public ArrayList<Particle>   P = new ArrayList<Particle>();

    /** List of Force objects. */
    public ArrayList<Force>      F = new ArrayList<Force>();

	/** Neighbor cache */
	public ArrayListParticle[] neighbors;
	public int numNeighbors = 0;

	/** Number of cells per dimension */
	public final double h = 0.1;
	public final int cellsPerLength = (int) Math.floor(dims / h);
//	public static final int numCells = cellsPerLength * cellsPerLength * cellsPerLength;


	/** Wrapper class so we can make an array of ArrayList<Particle>s
	 *  Because Java is a silly, silly language
	 *  */
	public class ArrayListParticle extends ArrayList<Particle> {}

	/** List of cells which keep track of particles */
	ArrayListParticle[][][] cells = new ArrayListParticle[cellsPerLength][cellsPerLength][cellsPerLength];

    /**
     * true iff prog has been initialized. This cannot be done in the
     * constructor because it requires a GL2 reference.
     */
    private boolean init = false;

    /** Filename of vertex shader source. */
    public static final String[] VERT_SOURCE = {"vert.glsl"};

    /** Filename of fragment shader source. */
    public static final String[] FRAG_SOURCE = {"frag.glsl"};

    /** The shader program used by the particles. */
    ShaderProgram prog;

	/** Cache variables */
	public static BlockingQueue<Vector3d> tempPool = new ArrayBlockingQueue<Vector3d>(numProcessors);

    /** Executor Service for concurrency */
    ExecutorService executor;


    /** Basic constructor. */
    public ParticleSystem() {
	    for (int i = 0; i < cellsPerLength; i++){
		    for (int j = 0; j < cellsPerLength; j++) {
			    for (int k = 0; k < cellsPerLength; k++) {
				    cells[i][j][k] = new ArrayListParticle();
			    }
		    }
	    }
	    neighbors = new ArrayListParticle[27];
	    Kernels.setH(h);
	    q = Kernels.poly6(0.03);

        for (int i = 0; i < numProcessors * 2; i++) {
            try {
                tempPool.put(new Vector3d());
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        executor = Executors.newFixedThreadPool(numProcessors);
        System.out.println("Available processors: " + numProcessors);
   }

    /**
     * Set up the GLSL program. This requires that the current directory (i.e. the package in which
     * this class resides) has a vertex and fragment shader.
     */
    public synchronized void init(GL2 gl) {
        if (init) return;

        prog = new ShaderProgram();
        ShaderCode vert_code = ShaderCode.create(gl, GL2ES2.GL_VERTEX_SHADER, 1, this.getClass(), VERT_SOURCE, false);
        ShaderCode frag_code = ShaderCode.create(gl, GL2ES2.GL_FRAGMENT_SHADER, 1, this.getClass(), FRAG_SOURCE, false);
        if (!prog.add(gl, vert_code, System.err) || !prog.add(gl, frag_code, System.err)) {
            System.err.println("WARNING: shader did not compile");
            prog.init(gl); // Initialize empty program
        } else {
            prog.link(gl, System.err);
        }

        init = true;
    }

    /** Adds a force object (until removed) */
    public synchronized void addForce(Force f) {
        F.add(f);
    }

    /** Useful for removing temporary forces, such as user-interaction
     * spring forces. */
    public synchronized void removeForce(Force f) {
        F.remove(f);
    }

    /** Creates particle and adds it to the particle system.
     * @param p0 Undeformed/material position.
     * @return Reference to new Particle.
     */
    public synchronized Particle createParticle(Point3d p0)
    {
        Particle newP = new Particle(p0);
        P.add(newP);
        return newP;
    }

    /**
     * Helper-function that computes nearest particle to the specified
     * (deformed) position.
     * @return Nearest particle, or null if no particles.
     */
    public synchronized Particle getNearestParticle(Point3d x)
    {
        Particle minP      = null;
        double   minDistSq = Double.MAX_VALUE;
        for(Particle particle : P) {
            double distSq = x.distanceSquared(particle.x);
            if(distSq < minDistSq) {
                minDistSq = distSq;
                minP = particle;
            }
        }
        return minP;
    }

    /** Moves all particles to undeformed/materials positions, and
     * sets all velocities to zero. Synchronized to avoid problems
     * with simultaneous calls to advanceTime(). */
    public synchronized void reset()
    {
        for(Particle p : P)  {
            p.x.set(p.x0);
            p.v.set(0,0,0);
            p.f.set(0,0,0);
            p.setHighlight(false);
        }
        time = 0;
    }

    /**
     * Incomplete/Debugging implementation of Forward-Euler
     * step. WARNING: Contains buggy debugging forces.
     */
    public synchronized void advanceTime(double dt)
    {

	    updateVandXstar(dt);
	    updateCells(false);
	    mainIteration(dt);
	    updateVandX(dt);

        time += dt;
    }

    /**
     * Displays Particle and Force objects.
     */
    public synchronized void display(GL2 gl)
    {
        for(Force force : F) {
            force.display(gl);
        }

        if(!init) init(gl);

        prog.useProgram(gl, true);

        for(Particle particle : P) {
            particle.display(gl);
        }

        prog.useProgram(gl, false);
    }

	public synchronized void updateCells(boolean useXStar) {
		int i, j, k;
		int prevCell[];
		for (Particle p : P){
			if (useXStar) {
				i = (int) Math.floor(p.xStar.x / h);
				j = (int) Math.floor(p.xStar.y / h);
				k = (int) Math.floor(p.xStar.z / h);
			} else {
				i = (int) Math.floor(p.x.x / h);
				j = (int) Math.floor(p.x.y / h);
				k = (int) Math.floor(p.x.z / h);
			}

			prevCell = p.cell;
			if (i != prevCell[0] || j!= prevCell[1] || k != prevCell[2]) {
				cells[prevCell[0]][prevCell[1]][prevCell[2]].remove(p);
				cells[i][j][k].add(p);
				prevCell[0] = i; prevCell[1] = j; prevCell[2] = k;
				p.cell = prevCell;
			}
		}
	}

	public synchronized void updateVandXstar(double dt){
		for (Particle p : P) {

            Vector3d temp1 = null;
            try {
                temp1 = tempPool.take();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            fExt(p.x, temp1);
			p.f.add(temp1);
			// v_i <= v_i + dt*fExt
			p.v.scaleAdd(dt, p.f, p.v);
			p.f.set(0,0,0);

			// x_i* <= x_i + dt*v_i
			p.xStar.scaleAdd(dt, p.v, p.x);

            try {
                tempPool.put(temp1);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
	}

	public synchronized void mainIteration(double dt) {

		for (int iter = 0; iter < maxIter; iter++) {
			// Calculate lambda for each particle

            Vector3d temp1 = null;
            try {
                temp1 = tempPool.take();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

			for (Particle p : P) {
				p.rho = 0;
				double delCSq = 0;
                temp1.set(0, 0, 0);

				// Calculate local density
				// Also calculate delC^2
				getNeighborCells(p);
				for (int i = 0; i < numNeighbors; i++) {
					for (Particle otherP : neighbors[i]) {
						otherP.rHat.sub(p.xStar, otherP.xStar);
						otherP.r = otherP.rHat.length();
						if (otherP.r > h || crossesBounds(p.x, otherP.x)) continue;

						p.rho += Kernels.poly6(otherP.r);

						// delCsq = sum_{i!=j}{abs(dW_ij)^2) + norm(sum_{i!=j}{dW_ij*rHat_ij}}^2
						// Kinda like a variance? But adding instead of subtracting?
						if (otherP.r != 0) {
							otherP.rHat.normalize();
							double dW = Kernels.delSpiky(otherP.r);
							temp1.scaleAdd(dW, otherP.rHat, temp1);
							delCSq += dW * dW;
						}
					}
				}
				double C = p.rho / rho_0 - 1;
				delCSq += temp1.lengthSquared();
				delCSq /= rho_0 * rho_0;
				delCSq += epsilon;
				p.lambda = -C / delCSq;
			}

			// Calculate deltaP for each particle
			for (Particle p : P) {
				getNeighborCells(p);
				p.dP.set(0, 0, 0);
				for (int i = 0; i < numNeighbors; i++) {
					for (Particle otherP : neighbors[i]) {
						otherP.rHat.set(p.xStar);
						otherP.rHat.sub(otherP.xStar);
						otherP.r = otherP.rHat.length();
						if (otherP.r > h || otherP.r == 0 || crossesBounds(p.x, otherP.x)) continue;
						otherP.rHat.normalize();
						double dW = Kernels.delSpiky(otherP.r);
						double s_corr = Kernels.poly6(otherP.r)/q;
						s_corr *= s_corr; s_corr *= s_corr;
						s_corr *= k;
						p.dP.scaleAdd((p.lambda + otherP.lambda + s_corr) * dW, otherP.rHat, p.dP);
					}
				}
				p.dP.scale(1 / rho_0);

				temp1.add(p.xStar, p.dP);

				if (temp1.x <= 0) {
					p.dP.x = -p.xStar.x;
					p.f.x = -p.v.x * p.m / dt;
				}
				if (temp1.y <= 0) {
					p.dP.y = -p.xStar.y;
					p.f.y = -p.v.y * p.m / dt;
				}
				if (temp1.z <= 0) {
					p.dP.z = -p.xStar.z;
					p.f.z = -p.v.z * p.m / dt;
				}
				if (temp1.x >= dims) {
					p.dP.x = dims - 0.001 - p.xStar.x;
					p.f.x = -p.v.x * p.m / dt;
				}
				if (temp1.y >= dims/2.0) {
					p.dP.y = dims/2.0 - 0.001 - p.xStar.y;
					p.f.y = -p.v.y * p.m / dt;
				}
				if (temp1.z >= dims/2.0) {
					p.dP.z = dims/2.0 - 0.001 - p.xStar.z;
					p.f.z = -p.v.z * p.m / dt;
				}

                if ((temp1.x > bounds1) && (p.x.x < bounds1) && ((p.x.z < dims * 0.45))) {
                    p.dP.x = bounds1 - 0.001 - p.xStar.x;
                    p.f.x = -p.v.x * p.m / dt;
                }


                if ((temp1.x < bounds1) && (p.x.x > bounds1) && ((p.x.z < dims * 0.45)) ) {
                    p.dP.x = bounds1 + 0.001 - p.xStar.x;
                    p.f.x = -p.v.x * p.m / dt;
                }

                if ((temp1.x > bounds1) && (p.x.x < bounds1) && ((p.x.z < dims * 0.45))) {
                    p.dP.x = bounds1 - 0.001 - p.xStar.x;
                    p.f.x = -p.v.x * p.m / dt;
                }

                if ((temp1.x < bounds1) && (p.x.x > bounds1) && ((p.x.z < dims * 0.45)) ) {
                    p.dP.x = bounds1 + 0.001 - p.xStar.x;
                    p.f.x = -p.v.x * p.m / dt;
                }

                if ((temp1.x > bounds2)) {
                    p.dP.x = bounds2 - 0.001 - p.xStar.x;
                    p.f.x = -p.v.x * p.m / dt;
                }

                try {
                    tempPool.put(temp1);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }

            }
			for (Particle p : P) {
				p.xStar.add(p.dP);
			}
		}

	}

	public synchronized void updateVandX(double dt) {
        for (Particle p : P) {
            p.v.sub(p.xStar, p.x);
            p.v.scale(1 / dt);
        }

        Vector3d temp1 = null;
        Vector3d temp2 = null;
        try {
            temp1 = tempPool.take();
            temp2 = tempPool.take();


            //XSPH Viscosity. Also calculate omegas
            for (Particle p : P) {
                getNeighborCells(p);
                temp1.set(0, 0, 0);
                p.omega.set(0, 0, 0);
                for (int i = 0; i < numNeighbors; i++) {
                    for (Particle otherP : neighbors[i]) {
                        otherP.rHat.sub(p.xStar, otherP.xStar);
                        otherP.r = otherP.rHat.length();
                        if (otherP.r < h && !crossesBounds(p.x, otherP.x)) {
                            temp2.sub(otherP.v, p.v);
                            temp1.scaleAdd(Kernels.poly6(otherP.r) / otherP.rho, temp2, temp1);
                            if (otherP.r != 0) {
                                otherP.rHat.normalize();
                                temp2.cross(temp2, otherP.rHat);
                                p.omega.scaleAdd(Kernels.delSpiky(otherP.r) / otherP.rho, temp2, p.omega);
                            }
                        }
                    }
                }

                p.vStar.scaleAdd(c, temp1, p.v);
            }

            // Vorticity confinement
            for (Particle p : P) {
                getNeighborCells(p);
                temp1.set(0, 0, 0);
                for (int i = 0; i < numNeighbors; i++) {
                    for (Particle otherP : neighbors[i]) {
                        otherP.rHat.sub(p.xStar, otherP.xStar);
                        otherP.r = otherP.rHat.length();
                        if (otherP.r < h && otherP.r != 0 && !crossesBounds(p.x, otherP.x)) {
                            otherP.rHat.normalize();
                            temp1.scaleAdd(otherP.omega.length() * Kernels.delSpiky(otherP.r) / otherP.rho, otherP.rHat, temp1);
                        }
                    }
                }
                if (temp1.length() != 0) temp1.normalize();
                temp1.cross(temp1, p.omega);
                p.f.scaleAdd(eVorticity, temp1, p.f);
            }

            for (Particle p : P) {
                p.v.set(p.vStar);
            }

            for (Particle p : P) {
                p.x.set(p.xStar);
            }

            tempPool.put(temp1);
            tempPool.put(temp2);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

	public synchronized int[] getCellIndices(Vector3d p) {
		int[] cellInd = new int[3];
		cellInd[0] = (int) Math.floor(p.x / h);
		cellInd[1] = (int) Math.floor(p.y / h);
		cellInd[2] = (int) Math.floor(p.z / h);
		return cellInd;
	}

	public synchronized void getNeighborCells(Particle p) {
		int[] cellInd = p.cell;
		numNeighbors = 0;
		for (int i = cellInd[0] - 1; i < cellInd[0] + 2; i++){
			for (int j = cellInd[1] - 1; j < cellInd[1] + 2; j++) {
				for (int k = cellInd[2] - 1; k < cellInd[2] + 2; k++) {
					if (i > -1 && j > -1 && k > -1 && i < cellsPerLength && j < cellsPerLength && k < cellsPerLength){
						neighbors[numNeighbors] = cells[i][j][k];
						numNeighbors++;
					}
				}
			}
		}
	}

	/** Not yet implemented. To be implemented if found necessary */
	public synchronized ArrayList<ArrayListParticle> getNeighborCells(Vector3d p) {
		int[] cellInd = getCellIndices(p);
		return null;
	}

	/** Position based external forces
	 *  This is actually acceleration, but all particles have uniform mass, so it's fine.
	 * */
	public void fExt(Point3d x, Vector3d force) {
		if (x.y < 0.1 && time >= 5) {
			double x0 = x.x - dims / 2.0;
			double z0 = x.z - dims / 2.0;
			double r = Math.sqrt(x0 * x0 + z0 * z0);
			double theta = Math.atan2(x0, z0);
			force.set(whirl * r * Math.cos(theta), -9.806, -whirl * r * Math.sin(theta));
		} else {
			force.set(0, -9.806, 0);
		}
    }

    public synchronized boolean crossesBounds(Tuple3d v1, Tuple3d v2){
        if ((v1.x - bounds1) * (v2.x - bounds1) <= 0 || (v1.x - bounds1) * (v2.x - bounds1) <= 0 ) return true;
        else return false;
    }

    public class Neighbor {
        int numNeighbors;

    }
}
