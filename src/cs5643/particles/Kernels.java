package cs5643.particles;

import sun.management.HotspotCompilationMBean;

import javax.vecmath.Vector3d;

/**
 * Created by Gagik on 2/24/2015.
 */
public class Kernels {
	private static double h = 0;
	private static double hCb = 0;
	private static double hSq = 0;
	private static double poly6Factor = 0;
	private static double spikyFactor = 0;
	private static double delSpikyFactor = 0;
	private Kernels(){}

	public static void setH(double H){
		h = H;
		hSq = h * h;
		hCb = h * h * h;
		poly6Factor = 315.0 / (64.0 * Math.PI * hCb * hCb * hCb);
		spikyFactor = 15.0 / (Math.PI * hCb * hCb);
		delSpikyFactor = -3 * spikyFactor;
	}

	public static double poly6(double r){
		if (r < h){
			double temp = hSq - r * r;
			return poly6Factor * temp * temp * temp;
		}
		else return 0;
	}

	public static double spiky(double r){
		if (r < h) {
			double temp = h - r;
			return spikyFactor * temp * temp * temp;
		} else return 0;
	}

	public static double delSpiky(double r){
		if (r < h) {
			double temp = h - r;
			return delSpikyFactor * temp * temp;
		} else return 0;
	}
}
