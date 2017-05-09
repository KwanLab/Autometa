package lu.uni.lcsb.vizbin.pca;

import no.uib.cipr.matrix.NotConvergedException;

public interface IPrincipleComponentAnalysis {
	void setup(int numSamples, int sampleSize);

	void addSample(double[] sampleData);

	void computeBasis(int numComponents) throws NotConvergedException;

	double[] sampleToEigenSpace(double[] sampleData);

}
