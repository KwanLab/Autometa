package lu.uni.lcsb.vizbin.graphics;

import java.awt.Polygon;

public class RegularPolygon extends Polygon {
	/**
	 * 
	 */
	private static final long	serialVersionUID	= 1L;

	public RegularPolygon(int x, int y, int r, int vertexCount) {
		this(x, y, r, vertexCount, 0);
	}

	public RegularPolygon(int x, int y, int r, int vertexCount, double startAngle) {
		super(getXCoordinates(x, y, r, vertexCount, startAngle), getYCoordinates(x, y, r, vertexCount, startAngle), vertexCount);
	}

	protected static int[] getXCoordinates(int x, int y, int r, int vertexCount, double startAngle) {
		int res[] = new int[vertexCount];
		double addAngle = 2 * Math.PI / vertexCount;
		double angle = startAngle;
		for (int i = 0; i < vertexCount; i++) {
			res[i] = (int) Math.round(r * Math.cos(angle)) + x;
			angle += addAngle;
		}
		return res;
	}

	protected static int[] getYCoordinates(int x, int y, int r, int vertexCount, double startAngle) {
		int res[] = new int[vertexCount];
		double addAngle = 2 * Math.PI / vertexCount;
		double angle = startAngle;
		for (int i = 0; i < vertexCount; i++) {
			res[i] = (int) Math.round(r * Math.sin(angle)) + y;
			angle += addAngle;
		}
		return res;
	}
}
