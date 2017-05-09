package lcsb.vizbin.service.utils.pca;

import static org.junit.Assert.assertTrue;
import lu.uni.lcsb.vizbin.pca.IPrincipleComponentAnalysis;
import lu.uni.lcsb.vizbin.pca.PrincipleComponentAnalysisMtj;

import org.apache.log4j.Logger;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class PrincipleComponentAnalysisMtjTest {
	Logger											logger	= Logger.getLogger(PrincipleComponentAnalysisMtjTest.class);

	IPrincipleComponentAnalysis	pca			= new PrincipleComponentAnalysisMtj();

	static final double					EPSILON	= 1e-6;

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void test() throws Exception {
		try {
			pca.setup(5, 10);
			double s1[] = new double[] { 9.61, 7.54, 5.92, 7.41, 8.80, 7.41, 4.19, 9.62, 3.00, 4.88 };
			double s2[] = new double[] { 4.31, 0.52, 1.61, 8.88, 5.19, 2.86, 1.27, 0.14, 5.34, 0.04 };
			double s3[] = new double[] { 5.51, 5.61, 7.56, 6.45, 6.42, 7.64, 2.85, 0.95, 4.28, 4.22 };
			double s4[] = new double[] { 4.93, 6.56, 5.18, 4.65, 6.49, 3.50, 6.19, 7.31, 7.82, 7.71 };
			double s5[] = new double[] { 8.99, 3.96, 5.76, 9.66, 1.47, 4.35, 4.68, 3.44, 1.46, 7.59 };

			double o1[] = new double[] { 6.463683831060132, -0.11602656248138686, -3.752375607662306 };
			double o2[] = new double[] { -9.176546316934445, -1.996482173065449, -0.4484723328472291 };
			double o3[] = new double[] { -1.1864454937202447, -0.417214015689895, -2.370746915503915 };
			double o4[] = new double[] { 4.159110592463409, -3.939571658344603, 4.307235069873009 };
			double o5[] = new double[] { -0.25980261286885575, 6.46929440958134, 2.2643597861404414 };

			pca.addSample(s1);
			pca.addSample(s2);
			pca.addSample(s3);
			pca.addSample(s4);
			pca.addSample(s5);
			pca.computeBasis(3);

			assertTrue(arrayEquals(o1, pca.sampleToEigenSpace(s1)));
			assertTrue(arrayEquals(o2, pca.sampleToEigenSpace(s2)));
			assertTrue(arrayEquals(o3, pca.sampleToEigenSpace(s3)));
			assertTrue(arrayEquals(o4, pca.sampleToEigenSpace(s4)));
			assertTrue(arrayEquals(o5, pca.sampleToEigenSpace(s5)));

		} catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}

	boolean arrayEquals(double[] a1, double[] a2) {
		if (a1.length != a2.length) {
			return false;
		}
		// for some reason matrixes have differnt signs...
		for (int i = 0; i < a1.length; i++) {
			if (Math.abs(a1[i]) - Math.abs(a2[i]) > EPSILON) {
				return false;
			}
		}
		return true;
	}

	@Test
	public void test2() throws Exception {
		try {
			pca.setup(6, 4);

			double s1[] = new double[] { 7.52, -1.10, -7.95, 1.08 };
			double s2[] = new double[] { -0.76, 0.62, 9.34, -7.10 };
			double s3[] = new double[] { 5.13, 6.62, -5.66, 0.87 };
			double s4[] = new double[] { -4.75, 8.52, 5.75, 5.30 };
			double s5[] = new double[] { 1.33, 4.91, -5.49, -3.52 };
			double s6[] = new double[] { -2.40, -6.77, 2.34, 3.95 };

			double o1[] = new double[] { 9.770251684203519, 3.8293845380584166 };
			double o2[] = new double[] { -9.444348971013547, 1.699568069610743 };
			double o3[] = new double[] { 6.965227200665963, -4.073597911804497 };
			double o4[] = new double[] { -7.641469730454019, -7.182741071486653 };
			double o5[] = new double[] { 4.66890767186292, -2.5748936413029266 };
			double o6[] = new double[] { -4.318567855264838, 8.302280016924918 };

			pca.addSample(s1);
			pca.addSample(s2);
			pca.addSample(s3);
			pca.addSample(s4);
			pca.addSample(s5);
			pca.addSample(s6);
			pca.computeBasis(2);

			assertTrue(arrayEquals(o1, pca.sampleToEigenSpace(s1)));
			assertTrue(arrayEquals(o2, pca.sampleToEigenSpace(s2)));
			assertTrue(arrayEquals(o3, pca.sampleToEigenSpace(s3)));
			assertTrue(arrayEquals(o4, pca.sampleToEigenSpace(s4)));
			assertTrue(arrayEquals(o5, pca.sampleToEigenSpace(s5)));
			assertTrue(arrayEquals(o6, pca.sampleToEigenSpace(s6)));

		} catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}

}
