package lu.uni.lcsb.vizbin.pca;

import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.factory.SingularValueDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;
import org.ejml.ops.SingularOps;

/**
 * <p>
 * The following is a simple example of how to perform basic principle component
 * analysis in EJML.
 * </p>
 * 
 * <p>
 * Principal Component Analysis (PCA) is typically used to develop a linear
 * model for a set of data (e.g. face images) which can then be used to test for
 * membership. PCA works by converting the set of data to a new basis that is a
 * subspace of the original set. The subspace is selected to maximize
 * information.
 * </p>
 * <p>
 * PCA is typically derived as an eigenvalue problem. However in this
 * implementation {@link org.ejml.factory.SingularValueDecomposition SVD} is
 * used instead because it will produce a more numerically stable solution.
 * Computation using EVD requires explicitly computing the variance of each
 * sample set. The variance is computed by squaring the residual, which can
 * cause loss of precision.
 * </p>
 * 
 * <p>
 * Usage:<br>
 * 1) call setup()<br>
 * 2) For each sample (e.g. an image ) call addSample()<br>
 * 3) After all the samples have been added call computeBasis()<br>
 * 4) Call sampleToEigenSpace() , eigenToSampleSpace() , errorMembership() ,
 * response()
 * </p>
 * 
 * @author Peter Abeles
 */
public class PrincipleComponentAnalysisEJML implements IPrincipleComponentAnalysis {

	// principle component subspace is stored in the rows
	private DenseMatrix64F	vt;

	// how many principle components are used
	private int							numComponents;

	// where the data is stored
	private DenseMatrix64F	a	= new DenseMatrix64F(1, 1);
	private int							sampleIndex;

	// mean values of each element across all the samples
	private double					mean[];

	// List<double[]> samples = new ArrayList<double[]>();

	public PrincipleComponentAnalysisEJML() {
	}

	/**
	 * Must be called before any other functions. Declares and sets up internal
	 * data structures.
	 * 
	 * @param numSamples
	 *          Number of samples that will be processed.
	 * @param sampleSize
	 *          Number of elements in each sample.
	 */
	public void setup(int numSamples, int sampleSize) {
		mean = new double[sampleSize];
		a.reshape(numSamples, sampleSize, false);
		sampleIndex = 0;
		numComponents = -1;
	}

	/**
	 * Adds a new sample of the raw data to internal data structure for later
	 * processing. All the samples must be added before computeBasis is called.
	 * 
	 * @param sampleData
	 *          Sample from original raw data.
	 */
	public void addSample(double[] sampleData) {
		if (a.getNumCols() != sampleData.length)
			throw new IllegalArgumentException("Unexpected sample size");
		if (sampleIndex >= a.getNumRows())
			throw new IllegalArgumentException("Too many samples");

		for (int i = 0; i < sampleData.length; i++) {
			a.set(sampleIndex, i, sampleData[i]);
		}
		// samples.add(sampleData);
		sampleIndex++;
	}

	/**
	 * Computes a basis (the principle components) from the most dominant
	 * eigenvectors.
	 * 
	 * @param numComponents
	 *          Number of vectors it will use to describe the data. Typically much
	 *          smaller than the number of elements in the input vector.
	 */
	public void computeBasis(int numComponents) {
		if (numComponents > a.getNumCols())
			throw new IllegalArgumentException("More components requested that the data's length.");
		if (sampleIndex != a.getNumRows())
			throw new IllegalArgumentException("Not all the data has been added");
		if (numComponents > sampleIndex)
			throw new IllegalArgumentException("More data needed to compute the desired number of components");

		this.numComponents = numComponents;

		// compute the mean of all the samples
		for (int i = 0; i < a.getNumRows(); i++) {
			for (int j = 0; j < mean.length; j++) {
				mean[j] += a.get(i, j);
			}
		}
		for (int j = 0; j < mean.length; j++) {
			mean[j] /= a.getNumRows();
		}

		// subtract the mean from the original data
		for (int i = 0; i < a.getNumRows(); i++) {
			for (int j = 0; j < mean.length; j++) {
				a.set(i, j, a.get(i, j) - mean[j]);
			}
		}

		// Compute SVD and save time by not computing U
		SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(a.numRows, a.numCols, false, true, false);
		if (!svd.decompose(a)) {
			throw new RuntimeException("SVD failed");
		}

		vt = svd.getV(null, true);
		DenseMatrix64F w = svd.getW(null);

		// Singular values are in an arbitrary order initially
		SingularOps.descendingOrder(null, false, w, vt, true);

		// strip off unneeded components and find the basis
		vt.reshape(numComponents, mean.length, true);
	}

	/**
	 * Returns a vector from the PCA's basis.
	 * 
	 * @param which
	 *          Which component's vector is to be returned.
	 * @return Vector from the PCA basis.
	 */
	public double[] getBasisVector(int which) {
		if (which < 0 || which >= numComponents) {
			throw new IllegalArgumentException("Invalid component");
		}

		DenseMatrix64F v = new DenseMatrix64F(1, a.numCols);
		CommonOps.extract(vt, which, which + 1, 0, a.numCols, v, 0, 0);

		return v.data;
	}

	/**
	 * Converts a vector from sample space into eigen space.
	 * 
	 * @param sampleData
	 *          Sample space data.
	 * @return Eigen space projection.
	 */
	public double[] sampleToEigenSpace(double[] sampleData) {
		if (sampleData.length != a.getNumCols())
			throw new IllegalArgumentException("Unexpected sample length");
		DenseMatrix64F mean = DenseMatrix64F.wrap(a.getNumCols(), 1, this.mean);

		DenseMatrix64F s = new DenseMatrix64F(a.getNumCols(), 1, true, sampleData);
		DenseMatrix64F r = new DenseMatrix64F(numComponents, 1);

		CommonOps.sub(s, mean, s);

		CommonOps.mult(vt, s, r);

		return r.data;
	}

	/**
	 * Converts a vector from eigen space into sample space.
	 * 
	 * @param eigenData
	 *          Eigen space data.
	 * @return Sample space projection.
	 */
	public double[] eigenToSampleSpace(double[] eigenData) {
		if (eigenData.length != numComponents)
			throw new IllegalArgumentException("Unexpected sample length");

		DenseMatrix64F s = new DenseMatrix64F(a.getNumCols(), 1);
		DenseMatrix64F r = DenseMatrix64F.wrap(numComponents, 1, eigenData);

		CommonOps.multTransA(vt, r, s);

		DenseMatrix64F mean = DenseMatrix64F.wrap(a.getNumCols(), 1, this.mean);
		CommonOps.add(s, mean, s);

		return s.data;
	}

	/**
	 * <p>
	 * The membership error for a sample. If the error is less than a threshold
	 * then it can be considered a member. The threshold's value depends on the
	 * data set.
	 * </p>
	 * <p>
	 * The error is computed by projecting the sample into eigenspace then
	 * projecting it back into sample space and
	 * </p>
	 * 
	 * @param sampleA
	 *          The sample whose membership status is being considered.
	 * @return Its membership error.
	 */
	public double errorMembership(double[] sampleA) {
		double[] eig = sampleToEigenSpace(sampleA);
		double[] reproj = eigenToSampleSpace(eig);

		double total = 0;
		for (int i = 0; i < reproj.length; i++) {
			double d = sampleA[i] - reproj[i];
			total += d * d;
		}

		return Math.sqrt(total);
	}

	/**
	 * Computes the dot product of each basis vector against the sample. Can be
	 * used as a measure for membership in the training sample set. High values
	 * correspond to a better fit.
	 * 
	 * @param sample
	 *          Sample of original data.
	 * @return Higher value indicates it is more likely to be a member of input
	 *         dataset.
	 */
	public double response(double[] sample) {
		if (sample.length != a.numCols)
			throw new IllegalArgumentException("Expected input vector to be in sample space");

		DenseMatrix64F dots = new DenseMatrix64F(numComponents, 1);
		DenseMatrix64F s = DenseMatrix64F.wrap(a.numCols, 1, sample);

		CommonOps.mult(vt, s, dots);

		return NormOps.normF(dots);
	}

	/*
	 * @Override public double[] sampleToEigenSpace(int sample) { return
	 * sampleToEigenSpace(samples.get(sample)); }
	 */
}