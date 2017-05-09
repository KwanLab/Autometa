package lu.uni.lcsb.vizbin.data;

import java.awt.geom.Point2D;
import java.util.HashMap;
import java.util.Map;

/**
 * Describes single sequence.
 * 
 * @author Piotr Gawron
 * 
 */
public class Sequence {

	/**
	 * Location of the sequence on the visualization layer (2d point).
	 */
	private Point2D									location;
	/**
	 * Identfier of the sequence.
	 */
	private String									id;

	/**
	 * Label identifier.
	 */
	private Integer									labelId		= 0;
	/**
	 * Label name.
	 */
	private String									labelName	= "";

	/**
	 * Dna sequence of the sequence.
	 */
	private String									dna				= "";

	/**
	 * Gc content.
	 */
	private Double									gc;

	/**
	 * Coverage.
	 */
	private Double									coverage;

	/**
	 * If <code>true</code> then this sequence is important (should be somehow
	 * marked on the visualization process).
	 */
	private Boolean									marker;

	/**
	 * Length of the sequence. Sometimes instead of the length the log(length) is
	 * stored here...
	 */
	private Double									length;

	/**
	 * Vector describing the sequence (normalized amount of kmers of different
	 * types). Key in this array is kmer string transformed into integer.
	 */
	private double[]								descVector;

	/**
	 * Description vector after clr. Key in this array is kmer string transformed
	 * into integer.
	 */
	private double[]								clrVector;
	/**
	 * Vector describing sequence after pca.
	 */
	private double[]								pcaVector;

	/**
	 * List of kmers amount (map value) for a given kmer length (map key).
	 */
	private Map<Integer, Integer[]>	kmersMap	= new HashMap<Integer, Integer[]>();

	/**
	 * @return the location
	 * @see #location
	 */
	public Point2D getLocation() {
		return location;
	}

	/**
	 * @param location
	 *          the location to set
	 * @see #location
	 */
	public void setLocation(Point2D location) {
		this.location = location;
	}

	/**
	 * @return the id
	 * @see #id
	 */
	public String getId() {
		return id;
	}

	/**
	 * @param id
	 *          the id to set
	 * @see #id
	 */
	public void setId(String id) {
		this.id = id;
	}

	/**
	 * @return the labelId
	 * @see #labelId
	 */
	public Integer getLabelId() {
		return labelId;
	}

	/**
	 * @param labelId
	 *          the labelId to set
	 * @see #labelId
	 */
	public void setLabelId(Integer labelId) {
		this.labelId = labelId;
	}

	/**
	 * @return the labelName
	 * @see #labelName
	 */
	public String getLabelName() {
		return labelName;
	}

	/**
	 * @param labelName
	 *          the labelName to set
	 * @see #labelName
	 */
	public void setLabelName(String labelName) {
		this.labelName = labelName;
	}

	/**
	 * @return the dna
	 * @see #dna
	 */
	public String getDna() {
		return dna;
	}

	/**
	 * @param dna
	 *          the dna to set
	 * @see #dna
	 */
	public void setDna(String dna) {
		this.dna = dna;
	}

	/**
	 * @return the gc
	 * @see #gc
	 */
	public Double getGc() {
		return gc;
	}

	/**
	 * @param gc
	 *          the gc to set
	 * @see #gc
	 */
	public void setGc(Double gc) {
		this.gc = gc;
	}

	/**
	 * @return the coverage
	 * @see #coverage
	 */
	public Double getCoverage() {
		return coverage;
	}

	/**
	 * @param coverage
	 *          the coverage to set
	 * @see #coverage
	 */
	public void setCoverage(Double coverage) {
		this.coverage = coverage;
	}

	/**
	 * @return the marker
	 * @see #marker
	 */
	public Boolean getMarker() {
		return marker;
	}

	/**
	 * @param marker
	 *          the marker to set
	 * @see #marker
	 */
	public void setMarker(Boolean marker) {
		this.marker = marker;
	}

	/**
	 * @return the length
	 * @see #length
	 */
	public Double getLength() {
		return length;
	}

	/**
	 * @param length
	 *          the length to set
	 * @see #length
	 */
	public void setLength(Double length) {
		this.length = length;
	}

	/**
	 * @return the descVector
	 * @see #descVector
	 */
	public double[] getDescVector() {
		return descVector;
	}

	/**
	 * @param descVector
	 *          the descVector to set
	 * @see #descVector
	 */
	public void setDescVector(double[] descVector) {
		this.descVector = descVector;
	}

	/**
	 * @return the clrVector
	 * @see #clrVector
	 */
	public double[] getClrVector() {
		return clrVector;
	}

	/**
	 * @param clrVector
	 *          the clrVector to set
	 * @see #clrVector
	 */
	public void setClrVector(double[] clrVector) {
		this.clrVector = clrVector;
	}

	/**
	 * @return the pcaVector
	 * @see #pcaVector
	 */
	public double[] getPcaVector() {
		return pcaVector;
	}

	/**
	 * @param pcaVector
	 *          the pcaVector to set
	 * @see #pcaVector
	 */
	public void setPcaVector(double[] pcaVector) {
		this.pcaVector = pcaVector;
	}

	/**
	 * @param kmerLength
	 *          length of the kmer
	 * @return the kmers for given kmer length
	 */
	public Integer[] getKmers(int kmerLength) {
		return kmersMap.get(kmerLength);
	}

	/**
	 * @param kmerLength
	 *          length of kmer
	 * @param list
	 *          list of kmers for given kmer length
	 */
	public void setKmers(int kmerLength, Integer[] list) {
		this.kmersMap.put(kmerLength, list);
	}

}
