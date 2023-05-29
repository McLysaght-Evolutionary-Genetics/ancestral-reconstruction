package analysis.macrosynteny;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.TreeMap;
import org.apache.commons.math3.special.Gamma;

public class MSynRec {
  static final boolean VERBOSE = Boolean.getBoolean("VERBOSE");
  static final boolean LOGSUMEXP = Boolean.getBoolean("LOGSUMEXP");

  static boolean ESTIMATE_ALPHA = Boolean.getBoolean("ESTIMATE_ALPHA");
  static boolean ESTIMATE_BETA = Boolean.getBoolean("ESTIMATE_BETA");
  static final int MIN_FQ_TO_AVOID_SUBOPTIMAL_INFERENCE = Integer.getInteger("MIN_FQ_TO_AVOID_SUBOPTIMAL_INFERENCE",
      -Integer.MAX_VALUE + 99999);
  static boolean FIXED_ALPHA = Boolean.getBoolean("FIXED_ALPHA");

  enum InitializationMethod {
    CLUSTERING, RANDOM_QX
  }

  InitializationMethod initializationMethod_ = InitializationMethod.CLUSTERING;

  enum InferenceMethod {
    CVB0, CVB, VBEM
  }

  InferenceMethod inferenceMethod_ = InferenceMethod.CVB0;

  public void setInitializationMethod(InitializationMethod initializationMethod) {
    initializationMethod_ = initializationMethod;
  }

  public void setInferenceMethod(InferenceMethod inferenceMethod) {
    inferenceMethod_ = inferenceMethod;
  }

  public void setEstimateHyperParamters(boolean alpha, boolean beta) {
    ESTIMATE_ALPHA = alpha;
    ESTIMATE_BETA = beta;
  }

  // for data io
  String[] postWgdSpeciesNameList_;

  // analysis constants
  int K_;
  // Number of postWGD species
  int T_;
  // Numbers of chromosomes in postWGD species
  int[] Ct_;
  // Maximum number of co-orthologs for each postWGD species
  int[] Dt_;

  // data
  ArrayList<Segment> segmentList_ = new ArrayList<Segment>();
  int S_;

  // algorithm parameters (convergence criteria, etc.)
  static final double THRESHOLD_ALPHA = 0.00001;
  static final double THRESHOLD_FQ = 0.00001;
  static final int VBEM_MIN_ITERATION = Integer.getInteger("VBEM_MIN_ITERATION", 1);
  static final int VBEM_MAX_ITERATION = Integer.getInteger("VBEM_MAX_ITERATION", 100);
  static final int VBEM_MAX_ITERATION_PER_SEGMENT = Integer.getInteger("VBEM_MAX_ITERATION_PER_SEGMENT", 100);
  static final int NEWTON_MAX_ITERATION = Integer.getInteger("NEWTON_MAX_ITERATION", 100);

  static final double CVB_CONVERGENCE_QX = 0.001;
  static final int CVB_MIN_ITERATION = Integer.getInteger("CVB_MIN_ITERATION", 1);
  static final int CVB_MAX_ITERATION = Integer.getInteger("CVB_MAX_ITERATION", 100);

  // qx initialization paramters
  static final double epsilon_ = 1e-10;

  // prior hyper-parameters
  // alpha_[k] := parameters for the prior K-dimentional dirichlet
  double[] alpha_;
  // beta_[t][c] := beta^{(k,t)}_c
  double[][] beta_;

  // variational parameters
  // vqx_[s][g][k] := \bar{\phi}_{s,g}(k)
  double[][][] vqx_;
  // valpha_[s][k] := variational parameters for the prior K-dimentinal dirichlet
  double[][] valpha_;
  // vbeta_[t][k][c] := variational \bar{\beta}^{(k,t)}_c. Note that the order of
  // indexing parameters is vb[t][k][c] not vb[k][t][c].
  double[][][] vbeta_;

  // digamma of variational parameters alpha and beta
  // dg_vask_[s][k] = digamma(valpha_sk)
  double[][] dg_vask_;
  // dg_sum_vask_[s] = digamma(sum_k valpha_sk)
  double[] dg_sum_vask_;
  // dg_vbtkc_[t][k][c] = digamma(vbeta_tkc)
  double[][][] dg_vbtkc_;
  // dg_sum_vbtkc_[t][k] = digamma(sum_c vbeta_tkc)
  double[][] dg_sum_vbtkc_;

  // logging
  String workingDirectoryPath_ = null;

  // information required for calculating BIC
  int numberOfOrthologs_ = 0;
  int numberOfGenes_ = 0;

  class Gene {
    String id_;
    // Y_[t][d] := Y_{s,g,t,d}
    int[][] Y_;

    public Gene(String line, int[] speciesIndexToRead) {
      // -1 means it doesn't remove empty columns at the end of each row
      String[] row = line.split("\t", -1);

      if (row.length < 1 + T_) {
        System.err.println("MSynRec.Gene(): Input file format error. T_=" + T_ + "\nline=" + line);
      }

      id_ = row[0];
      Y_ = new int[T_][];

      for (int t = 0; t < T_; ++t) {
        int Dt = Dt_[t];
        int Ct = Ct_[t];
        Y_[t] = new int[Dt];
        String orthologList = row[1 + speciesIndexToRead[t]];

        // no orthologs
        if (orthologList.length() == 0) {
          continue;
        }

        String[] orthologs = orthologList.split(",");

        for (int j = 0; j < orthologs.length && j < Dt; ++j) {
          int c = Integer.parseInt(orthologs[j]);

          if (c >= Ct) {
            System.err.println(
                "MSynRec.Gene(): Input file format error. Ct[" + (t) + "]=" + Ct + "<=" + c + "\nline=" + line);
            System.err.println(t + "\t" + speciesIndexToRead[t]);
          }

          if (c > 0) {
            ++numberOfOrthologs_;
          }

          // =0: no co-ortholog, =1,...,Ct: chromosome numbers,
          Y_[t][j] = 1 + c;
        }
      }

      ++numberOfGenes_;
    }

    /**
     * Returns the chromosome index of the d-th co-ortholog in teleost t.
     * 
     * @param t
     * @param d
     * @return
     *         -1: no co-ortholog
     *         0,...,Ct-1: chromosome number of the ortholog
     */
    public int getYsg(int t, int d) {
      return Y_[t][d] - 1;
    }

    public String getID() {
      return id_;
    }
  }

  class Segment implements Comparable<Segment> {
    String id_;
    ArrayList<Gene> geneList_;

    public Segment(String id, ArrayList<Gene> geneList) {
      id_ = id;
      geneList_ = geneList;
    }

    public Gene getGene(int i) {
      return geneList_.get(i);
    }

    public int getLength() {
      return geneList_.size();
    }

    public String getID() {
      return id_;
    }

    public int[][] getOrthologCounts() {
      // initialize
      int[][] counts = new int[T_][];

      for (int t = 0; t < T_; ++t) {
        counts[t] = new int[Ct_[t]];
      }

      // count
      for (Gene g : geneList_) {
        for (int t = 0; t < T_; ++t) {
          for (int d = 0; d < Dt_[t]; ++d) {
            int chr = g.getYsg(t, d);

            // no co-ortholog
            if (chr == -1) {
              continue;
            }

            ++counts[t][chr];
          }
        }
      }

      return counts;
    }

    public int getOrthologCountTotal() {
      int[][] c = getOrthologCounts();
      int total = 0;

      for (int i = 0; i < c.length; ++i) {
        int[] d = c[i];

        for (int j = 0; j < d.length; ++j) {
          total += d[j];
        }
      }

      return total;
    }

    public int compareTo(Segment arg0) {
      return this.getOrthologCountTotal() - arg0.getOrthologCountTotal();
    }
  }

  public MSynRec(String[] postWgdSpeciesNameList, int K, int T, int[] Ct, int[] Dt) {
    postWgdSpeciesNameList_ = postWgdSpeciesNameList;
    K_ = K;
    T_ = T;
    Ct_ = Ct;
    Dt_ = Dt;
  }

  public void setWorkingDirectory(String path) {
    this.workingDirectoryPath_ = path;
  }

  public void setPrior(double[] alpha, double[][] beta) {
    alpha_ = alpha;
    beta_ = beta;
  }

  public void setPrior(double alpha, double beta) {
    alpha_ = getInitialAlphaPrior(alpha);
    beta_ = getInitialBetaPrior(beta);
  }

  private double[] getInitialAlphaPrior(double symmetricAlpha) {
    double[] prior = new double[K_];

    for (int k = 0; k < K_; ++k) {
      prior[k] = symmetricAlpha;
    }

    return prior;
  }

  private double[][] getInitialBetaPrior(double symmetricBeta) {
    double[][] prior = new double[T_][];

    for (int t = 0; t < T_; ++t) {
      prior[t] = new double[Ct_[t]];

      for (int c = 0; c < Ct_[t]; ++c) {
        prior[t][c] = symmetricBeta;
      }
    }

    return prior;
  }

  private void initializeVariationalAlpha(int S_) {
    valpha_ = new double[S_][];
    dg_vask_ = new double[S_][];
    dg_sum_vask_ = new double[S_];

    for (int s = 0; s < S_; ++s) {
      valpha_[s] = new double[K_];
      dg_vask_[s] = new double[K_];
    }
  }

  private void initializeVariationalAlphaForSegment(int seg, double[] Alpha) {
    dg_sum_vask_[seg] = 0;

    for (int k = 0; k < K_; ++k) {
      valpha_[seg][k] = Alpha[k];
      dg_vask_[seg][k] = MyGamma_digamma(valpha_[seg][k]);
      dg_sum_vask_[seg] += valpha_[seg][k];
    }

    dg_sum_vask_[seg] = MyGamma_digamma(dg_sum_vask_[seg]);
  }

  private void initializeVariationalBeta() {
    // Random r = new Random();
    vbeta_ = new double[T_][][];
    dg_vbtkc_ = new double[T_][][];
    dg_sum_vbtkc_ = new double[T_][];

    for (int t = 0; t < T_; ++t) {
      vbeta_[t] = new double[K_][];
      dg_vbtkc_[t] = new double[K_][];
      dg_sum_vbtkc_[t] = new double[K_];

      for (int k = 0; k < K_; ++k) {
        vbeta_[t][k] = new double[Ct_[t]];
        dg_vbtkc_[t][k] = new double[Ct_[t]];

        for (int c = 0; c < Ct_[t]; ++c) {
          // vbeta_[t][k][c] = symmetricBeta;
          // vbeta_[t][k][c] = 0.5+r.nextDouble();//Uniform[0.5,1.5]
          vbeta_[t][k][c] = this.beta_[t][c];
          dg_vbtkc_[t][k][c] = MyGamma_digamma(vbeta_[t][k][c]);
          dg_sum_vbtkc_[t][k] += vbeta_[t][k][c];
        }

        dg_sum_vbtkc_[t][k] = MyGamma_digamma(dg_sum_vbtkc_[t][k]);
      }
    }
  }

  private void initializeVariationalQX(ArrayList<Segment> segmentList_) {
    S_ = segmentList_.size();
    Random r = new Random();
    vqx_ = new double[S_][][];

    for (int s = 0; s < S_; ++s) {
      Segment seg = segmentList_.get(s);
      int Gs = seg.getLength();
      vqx_[s] = new double[Gs][];

      for (int g = 0; g < Gs; ++g) {
        vqx_[s][g] = new double[K_];
        double sum = 0;

        for (int k = 0; k < K_; ++k) {
          vqx_[s][g][k] = 1;

          if (initializationMethod_ == InitializationMethod.RANDOM_QX) {
            vqx_[s][g][k] = r.nextDouble();
          }

          sum += vqx_[s][g][k];
        }

        for (int k = 0; k < K_; ++k) {
          vqx_[s][g][k] /= sum;

          if (LOGSUMEXP) {
            vqx_[s][g][k] = Math.log(vqx_[s][g][k]);
          }
        }
      }
    }
  }

  private void initializeByClustering() throws IOException {
    TreeMap<String, Integer> seg2cluster = null;

    InitialSegmentClustering isc = new InitialSegmentClustering(this.segmentList_, T_, Ct_);
    seg2cluster = isc.computeClusters(K_, workingDirectoryPath_);

    for (int s = 0; s < S_; ++s) {
      Segment seg = this.segmentList_.get(s);
      int Gs = seg.getLength();
      Integer clusterIndex = seg2cluster.get(seg.getID());

      // initialize vqx_
      // small segments that are excluded from initialization are initialized by the
      // uniform distribution
      if (clusterIndex == null) {
        for (int g = 0; g < Gs; ++g) {
          for (int k = 0; k < K_; ++k) {
            vqx_[s][g][k] = 1 / (double) K_;
          }
        }
      } else {
        for (int g = 0; g < Gs; ++g) {
          for (int k = 0; k < K_; ++k) {
            if (clusterIndex.intValue() == k) {
              vqx_[s][g][k] = 1 - (K_ - 1) * epsilon_;
            } else
              vqx_[s][g][k] = epsilon_;
          }
        }
      }
      if (LOGSUMEXP) {
        for (int g = 0; g < Gs; ++g) {
          for (int k = 0; k < K_; ++k) {
            vqx_[s][g][k] = Math.log(vqx_[s][g][k]);
          }
        }
      }
    }
  }

  /**
   * Reads the parameter values S, G_t, and Y_{s,g,t,d} from the input file.
   * T and D_t are assumed to be given.
   * file format:
   * segID1 \n
   * geneID1 \t ortholog1,ortholog2 \t ortholog1 \t ...
   * geneID2 \t ...
   * (blank line)
   * segID2 \n
   * ...
   */
  public void readSyntenyData(String filename) throws IOException {
    readSyntenyFileAndMakeSegmentData(filename, this.segmentList_);
    this.S_ = this.segmentList_.size();
  }

  private void readSyntenyFileAndMakeSegmentData(String filename, ArrayList<Segment> segmentList) throws IOException {
    BufferedReader br = new BufferedReader(new FileReader(filename));
    String line;
    ArrayList<Gene> geneList = new ArrayList<Gene>();
    String segID = null;
    String speciesListInfo = new String("#species=");
    int[] postWgdSpeciesIndex = new int[T_];

    for (int t = 0; t < T_; ++t) {
      postWgdSpeciesIndex[t] = t;
    }

    while ((line = br.readLine()) != null) {
      // comment line
      if (line.isEmpty() == false && line.charAt(0) == '#') {
        if (line.startsWith(speciesListInfo)) {
          String speciesInfo = line.substring(speciesListInfo.length());
          String[] speciesList = speciesInfo.split(",");
          HashMap<String, Integer> species2index = new HashMap<String, Integer>();

          for (int s = 0; s < speciesList.length; ++s) {
            species2index.put(speciesList[s], Integer.valueOf(s));
          }

          for (int i = 0; i < T_; ++i) {
            Integer speciesIndexInInputFile = species2index.get(postWgdSpeciesNameList_[i]);

            if (speciesIndexInInputFile == null) {
              System.err.println("MSynRec.readSyntenyData(): Post-WGD species (" + postWgdSpeciesNameList_[i]
                  + ") not found in the input file. " + filename);
              System.exit(1);
            }

            postWgdSpeciesIndex[i] = speciesIndexInInputFile.intValue();
          }
        }

        continue;
      } else if (segID == null) {
        segID = line;
        continue;
      }

      if (line.length() == 0) {
        // got a complete segment
        // Skip if the segment has no genes.
        if (geneList.isEmpty() == false) {
          Segment seg = new Segment(segID, geneList);
          segmentList.add(seg);

          geneList = new ArrayList<Gene>();
        }

        segID = null;

        continue;
      }

      Gene gene = new Gene(line, postWgdSpeciesIndex);
      geneList.add(gene);
    }

    br.close();

    // Add the last segment
    // Skip if the segment has no genes.
    if (geneList.isEmpty() == false) {
      Segment seg = new Segment(segID, geneList);
      segmentList.add(seg);
    }
  }

  private double computeFq() {
    double FXq = 0;
    double FXUq = 0;
    double FXVq = 0;

    double sum_alpha_k = 0;

    if (alpha_ == null) {
      System.err.println("ERROR in computeFq(): Hyper-parameters not initialized? Please call setPrior(alpha,beta)");
    }

    for (int k = 0; k < K_; ++k) {
      FXUq -= MyGamma_logGamma(alpha_[k]);
      sum_alpha_k += alpha_[k];
    }

    FXUq += MyGamma_logGamma(sum_alpha_k);
    FXUq *= S_;

    for (int s = 0; s < S_; ++s) {
      Segment seg = this.segmentList_.get(s);
      int Gs = seg.getLength();

      for (int g = 0; g < Gs; ++g) {
        for (int k = 0; k < K_; ++k) {
          if (LOGSUMEXP) {
            FXq += Math.exp(vqx_[s][g][k]) * vqx_[s][g][k];
          } else {
            FXq += vqx_[s][g][k] * Math.log(vqx_[s][g][k]);
          }
        }
      }

      for (int k = 0; k < K_; ++k) {
        FXUq += MyGamma_logGamma(valpha_[s][k]);
      }

      FXUq -= MyGamma_logGamma(Gs + sum_alpha_k);
    }

    // compute FXVq
    for (int t = 0; t < T_; ++t) {
      for (int k = 0; k < K_; ++k) {
        double sum_vbeta_ktc = 0;

        for (int c = 0; c < Ct_[t]; ++c) {
          FXVq += MyGamma_logGamma(vbeta_[t][k][c]);
          sum_vbeta_ktc += vbeta_[t][k][c];
        }

        FXVq -= MyGamma_logGamma(sum_vbeta_ktc);
      }

      double sum_beta_tc = 0;

      for (int c = 0; c < Ct_[t]; ++c) {
        FXVq -= K_ * MyGamma_logGamma(beta_[t][c]);
        sum_beta_tc += beta_[t][c];
      }

      FXVq += K_ * MyGamma_logGamma(sum_beta_tc);
    }

    double Fq = -FXq + FXUq + FXVq;

    return Fq;
  }

  private double get_vqx_sgk(int s, int g, int k) {
    double value = 0;

    if (LOGSUMEXP) {
      value = Math.exp(vqx_[s][g][k]);
    } else {
      value = vqx_[s][g][k];
    }

    return value;
  }

  private double[] computeEqNsk(int s) {
    Segment seg = this.segmentList_.get(s);
    double[] EqNsk = new double[K_];

    for (int g = 0; g < seg.getLength(); ++g) {
      for (int k = 0; k < K_; ++k) {
        EqNsk[k] += get_vqx_sgk(s, g, k);
      }
    }

    return EqNsk;
  }

  interface AbstractELogGammaSum {
    public void add(int coeff, double p);

    public double getValue(double a);
  }

  interface ELogSumBernoulli extends AbstractELogGammaSum {
    public void add(int coeff, double p);

    public void remove(int coeff, double p);

    public double getValue(double a);
  }

  class ELogSumBernoulliTaylorApproximation1stOrder implements ELogSumBernoulli {
    double ex_ = 0;

    public ELogSumBernoulliTaylorApproximation1stOrder() {
    }

    public void add(int coeff, double p) {
      ex_ += coeff * p;
    }

    public void remove(int coeff, double p) {
      ex_ -= coeff * p;

      if (ex_ < 0) {
        ex_ = 0;
      }
    }

    public double getValue(double a) {
      if (ex_ > 0) {
        return Math.log(a + ex_);
      } else {
        return Math.log(a);
      }
    }
  }

  class ELogSumBernoulliTaylorApproximation2ndOrder implements ELogSumBernoulli {
    double ex_;
    double varx_;

    public ELogSumBernoulliTaylorApproximation2ndOrder() {
    }

    public void add(int coeff, double p) {
      ex_ += coeff * p;
      varx_ += coeff * coeff * p * (1 - p);
    }

    public void remove(int coeff, double p) {
      ex_ -= coeff * p;
      varx_ -= coeff * coeff * p * (1 - p);

      if (ex_ < 0) {
        ex_ = 0;
      }

      if (varx_ < 0) {
        varx_ = 0;
      }
    }

    public double getValue(double a) {
      double value = Math.log(a + ex_);
      if (ex_ > 0) {
        value -= varx_ / 2.0 / ex_ / ex_;
      }

      if (value < Math.log(epsilon_)) {
        value = Math.log(epsilon_);
      }

      return value;
    }
  }

  public double computeCVB() throws IOException {
    this.initializeVariationalQX(this.segmentList_);

    switch (initializationMethod_) {
      case CLUSTERING:
        initializeByClustering();
        break;
      case RANDOM_QX:
        throw new UnsupportedOperationException("unimplemented");
    }

    // initialize expectations
    ELogSumBernoulli[][][] ENktc = new ELogSumBernoulli[K_][T_][];
    ELogSumBernoulli[][] ENkt = new ELogSumBernoulli[K_][T_];

    for (int k = 0; k < K_; ++k) {
      ENktc[k] = new ELogSumBernoulli[T_][];
      ENkt[k] = new ELogSumBernoulli[T_];

      for (int t = 0; t < T_; ++t) {
        ENktc[k][t] = new ELogSumBernoulli[Ct_[t]];

        for (int c = 0; c < Ct_[t]; ++c) {
          // Choose E[log(a+\sum_i X_i)] approx method
          switch (inferenceMethod_) {
            case CVB0:
              ENktc[k][t][c] = new ELogSumBernoulliTaylorApproximation1stOrder();
              break;
            case CVB:
              ENktc[k][t][c] = new ELogSumBernoulliTaylorApproximation2ndOrder();
              break;
            case VBEM:
              throw new UnsupportedOperationException("unimplemented");
          }
        }

        // Choose E[log(a+\sum_i X_i)] approx method
        switch (inferenceMethod_) {
          case CVB0:
            ENkt[k][t] = new ELogSumBernoulliTaylorApproximation1stOrder();
            break;
          case CVB:
            ENkt[k][t] = new ELogSumBernoulliTaylorApproximation2ndOrder();
            break;
          case VBEM:
            throw new UnsupportedOperationException("unimplemented");
        }
      }
    }

    // compute expectation of the number of orthologs in chr=c and teleost=t that
    // are orthologous to genes with anc=k.
    for (int s = 0; s < S_; ++s) {
      Segment seg = this.segmentList_.get(s);

      for (int g = 0; g < seg.getLength(); ++g) {
        Gene gene = seg.getGene(g);

        for (int t = 0; t < T_; ++t) {
          int[] orthologCount_sgtc = new int[Ct_[t]];
          int orthologCount_sgt = 0;

          for (int d = 0; d < Dt_[t]; ++d) {
            int c = gene.getYsg(t, d);
            if (c == -1) {
              break;
            }

            ++orthologCount_sgtc[c];
            ++orthologCount_sgt;
          }

          for (int k = 0; k < K_; ++k) {
            double vqxsgk = get_vqx_sgk(s, g, k);

            for (int c = 0; c < Ct_[t]; ++c) {
              if (orthologCount_sgtc[c] == 0) {
                continue;
              }

              ENktc[k][t][c].add(orthologCount_sgtc[c], vqxsgk);
            }

            if (orthologCount_sgt == 0) {
              continue;
            }

            ENkt[k][t].add(orthologCount_sgt, vqxsgk);
          }
        }
      }
    }

    // update algorithm
    int iteration = 0;
    boolean convergence = false;
    double prevFq = -Double.MAX_VALUE;
    double smoothingFactor = 0.8;

    while (convergence == false) {
      // CVB step
      double diffqX = 0;
      double alpha_additional = 0;

      if (FIXED_ALPHA == false && iteration < CVB_MAX_ITERATION - 3) {
        alpha_additional = Math.pow(smoothingFactor, iteration);
      }

      // if(iteration>0){
      for (int s = 0; s < S_; ++s) {
        Segment seg = this.segmentList_.get(s);
        // For each segment s, compute expectation of the number of genes with anc=k in
        // segment s.
        ELogSumBernoulli[] ENsk = new ELogSumBernoulli[K_];

        for (int k = 0; k < K_; ++k) {
          // Choose E[log(a+\sum_i X_i)] approx method
          switch (inferenceMethod_) {
            case CVB0:
              ENsk[k] = new ELogSumBernoulliTaylorApproximation1stOrder();
              break;
            case CVB:
              ENsk[k] = new ELogSumBernoulliTaylorApproximation2ndOrder();
              break;
            case VBEM:
              throw new UnsupportedOperationException("unimplemented");
          }

          for (int g = 0; g < seg.getLength(); ++g) {
            double vqxsgk = get_vqx_sgk(s, g, k);
            ENsk[k].add(1, vqxsgk);
          }
        }

        for (int g = 0; g < seg.getLength(); ++g) {
          Gene gene = seg.getGene(g);
          double proportionalConstant = 0;
          double logproportionalConstant = Double.NaN;
          double[] prev_vqxsgk = new double[K_];

          for (int k = 0; k < K_; ++k) {
            // Save prev values for calculating diff=\sum|old-new| after updating.
            prev_vqxsgk[k] = vqx_[s][g][k];
            double vqxsgk = get_vqx_sgk(s, g, k);
            ENsk[k].remove(1, vqxsgk);
            double new_logqxsgk = ENsk[k].getValue(alpha_[k] + alpha_additional);

            for (int t = 0; t < T_; ++t) {
              int[] orthologCount_sgtc = new int[Ct_[t]];
              int orthologCount_sgt = 0;

              for (int d = 0; d < Dt_[t]; ++d) {
                int c = gene.getYsg(t, d);

                if (c == -1) {
                  break;
                }

                ++orthologCount_sgtc[c];
                ++orthologCount_sgt;
              }

              for (int c = 0; c < Ct_[t]; ++c) {
                if (orthologCount_sgtc[c] == 0) {
                  continue;
                }

                ENktc[k][t][c].remove(orthologCount_sgtc[c], vqxsgk);
              }
              if (orthologCount_sgt == 0) {
                continue;
              }

              ENkt[k][t].remove(orthologCount_sgt, vqxsgk);

              double sumBeta = 0;

              for (int c = 0; c < Ct_[t]; ++c) {
                sumBeta += beta_[t][c];

                for (int m = 0; m < orthologCount_sgtc[c]; ++m) {
                  new_logqxsgk += ENktc[k][t][c].getValue(beta_[t][c] + m + Math.pow(smoothingFactor, iteration));
                }
              }

              for (int m = 0; m < orthologCount_sgt; ++m) {
                new_logqxsgk -= ENkt[k][t].getValue(sumBeta + m + Math.pow(smoothingFactor, iteration));
              }
            }

            // save the new vqx_[s][g][k] value.
            if (LOGSUMEXP) {
              vqx_[s][g][k] = new_logqxsgk;
              logproportionalConstant = logsumexp(logproportionalConstant, vqx_[s][g][k]);
            } else {
              vqx_[s][g][k] = Math.exp(new_logqxsgk);
              proportionalConstant += vqx_[s][g][k];
            }
          }

          if (LOGSUMEXP) {
            for (int k = 0; k < K_; ++k) {
              vqx_[s][g][k] -= logproportionalConstant;
            }
          } else {
            for (int k = 0; k < K_; ++k) {
              vqx_[s][g][k] /= proportionalConstant;
            }
          }

          // Restore expectations and calculate diff
          for (int k = 0; k < K_; ++k) {
            double new_vqxsgk = get_vqx_sgk(s, g, k);
            ENsk[k].add(1, new_vqxsgk);

            for (int t = 0; t < T_; ++t) {
              int[] orthologCount_sgtc = new int[Ct_[t]];
              int orthologCount_sgt = 0;

              for (int d = 0; d < Dt_[t]; ++d) {
                int c = gene.getYsg(t, d);
                if (c == -1) {
                  break;
                }

                ++orthologCount_sgtc[c];
                ++orthologCount_sgt;
              }

              for (int c = 0; c < Ct_[t]; ++c) {
                if (orthologCount_sgtc[c] == 0) {
                  continue;
                }

                ENktc[k][t][c].add(orthologCount_sgtc[c], new_vqxsgk);
              }

              if (orthologCount_sgt == 0) {
                continue;
              }

              ENkt[k][t].add(orthologCount_sgt, new_vqxsgk);
            }

            if (LOGSUMEXP) {
              diffqX += Math.abs(Math.exp(prev_vqxsgk[k]) - new_vqxsgk);
            } else {
              diffqX += Math.abs(prev_vqxsgk[k] - new_vqxsgk);
            }
          }

        }
      }

      ++iteration;

      if (CVB_MIN_ITERATION <= iteration && diffqX < CVB_CONVERGENCE_QX) {
        convergence = true;
      }

      if (iteration >= CVB_MAX_ITERATION) {
        convergence = true;
      }
    }

    this.initializeVariationalAlpha(this.S_);

    for (int s = 0; s < S_; ++s) {
      // For each segment s, compute expectation of the number of genes with anc=k in
      // segment s.
      double[] EqNsk = computeEqNsk(s);

      for (int k = 0; k < K_; ++k) {
        valpha_[s][k] = alpha_[k] + EqNsk[k];
      }
    }

    this.initializeVariationalBeta();

    for (int t = 0; t < T_; ++t) {
      for (int k = 0; k < K_; ++k) {
        for (int c = 0; c < Ct_[t]; ++c) {
          vbeta_[t][k][c] = beta_[t][c];

          for (int s = 0; s < S_; ++s) {
            Segment seg = this.segmentList_.get(s);

            for (int g = 0; g < seg.getLength(); ++g) {
              Gene gene = seg.getGene(g);
              double vqxsgk = get_vqx_sgk(s, g, k);

              for (int d = 0; d < Dt_[t]; ++d) {
                int cc = gene.getYsg(t, d);

                if (cc == -1) {
                  break;
                }

                if (cc == c) {
                  vbeta_[t][k][c] += vqxsgk;
                }
              }
            }
          }

          dg_vbtkc_[t][k][c] = MyGamma_digamma(vbeta_[t][k][c]);
          dg_sum_vbtkc_[t][k] += vbeta_[t][k][c];
        }

        dg_sum_vbtkc_[t][k] = MyGamma_digamma(dg_sum_vbtkc_[t][k]);
      }
    }

    prevFq = computeFq();

    return prevFq;
  }

  public double computeVBEM() throws IOException {
    // intialize variational parameters
    this.initializeVariationalBeta();
    this.initializeVariationalAlpha(this.S_);
    this.initializeVariationalQX(this.segmentList_);

    switch (initializationMethod_) {
      case CLUSTERING:
        initializeByClustering();
        break;
      case RANDOM_QX:
        throw new UnsupportedOperationException("unimplemented");
    }

    // update algorithm
    int iteration = 0;
    boolean convergence = false;
    double prevFq = -Double.MAX_VALUE;

    while (convergence == false) {
      // VB-E step
      if (this.initializationMethod_ != InitializationMethod.CLUSTERING || iteration > 0) {
        for (int s = 0; s < S_; ++s) {
          initializeVariationalAlphaForSegment(s, this.alpha_);
          int iteration_segment = 0;
          double diffAlphaToCheckConvergence = 0;
          Segment seg = this.segmentList_.get(s);

          do {
            // VB-E step for vqx_
            for (int g = 0; g < seg.getLength(); ++g) {
              Gene gene = seg.getGene(g);
              double proportionalConstant = 0;
              double logproportionalConstant = Double.NaN;

              for (int k = 0; k < K_; ++k) {
                double logqxsgk = this.dg_vask_[s][k];

                for (int t = 0; t < T_; ++t) {
                  for (int d = 0; d < Dt_[t]; ++d) {
                    if (gene.getYsg(t, d) == -1) {
                      break;
                    }

                    logqxsgk += this.dg_vbtkc_[t][k][gene.getYsg(t, d)];
                    logqxsgk -= this.dg_sum_vbtkc_[t][k];
                  }
                }

                if (LOGSUMEXP) {
                  vqx_[s][g][k] = logqxsgk;
                  logproportionalConstant = logsumexp(logproportionalConstant, vqx_[s][g][k]);
                } else {
                  vqx_[s][g][k] = Math.exp(logqxsgk);
                  proportionalConstant += vqx_[s][g][k];
                }
              }

              if (LOGSUMEXP) {
                for (int k = 0; k < K_; ++k) {
                  vqx_[s][g][k] -= logproportionalConstant;
                }
              } else {
                for (int k = 0; k < K_; ++k) {
                  vqx_[s][g][k] /= proportionalConstant;
                }
              }
            }

            // VB-M step for valpha_
            this.dg_sum_vask_[s] = 0;

            for (int k = 0; k < K_; ++k) {
              double oldAlpha = valpha_[s][k];
              valpha_[s][k] = alpha_[k];

              for (int g = 0; g < seg.getLength(); ++g) {
                if (LOGSUMEXP) {
                  valpha_[s][k] += Math.exp(vqx_[s][g][k]);
                } else {
                  valpha_[s][k] += vqx_[s][g][k];
                }
              }

              // update digamma(valpha)
              this.dg_vask_[s][k] = MyGamma_digamma(valpha_[s][k]);
              // dg_sum_vask_[s] stores summand for every k
              this.dg_sum_vask_[s] += valpha_[s][k];

              diffAlphaToCheckConvergence += Math.abs(valpha_[s][k] - oldAlpha) / (double) K_;
            }

            this.dg_sum_vask_[s] = MyGamma_digamma(this.dg_sum_vask_[s]);
          } while (diffAlphaToCheckConvergence >= THRESHOLD_ALPHA
              && ++iteration_segment < VBEM_MAX_ITERATION_PER_SEGMENT);
        }
      }

      // VB-M step for vbeta_
      for (int k = 0; k < K_; ++k) {
        for (int t = 0; t < T_; ++t) {
          this.dg_sum_vbtkc_[t][k] = 0;

          for (int c = 0; c < Ct_[t]; ++c) {
            vbeta_[t][k][c] = beta_[t][c];

            for (int s = 0; s < S_; ++s) {
              Segment seg = this.segmentList_.get(s);

              for (int g = 0; g < seg.getLength(); ++g) {
                Gene gene = seg.getGene(g);

                for (int d = 0; d < Dt_[t]; ++d) {
                  if (gene.getYsg(t, d) != c) {
                    continue;
                  }
                  if (LOGSUMEXP) {
                    vbeta_[t][k][c] += Math.exp(vqx_[s][g][k]);
                  } else {
                    vbeta_[t][k][c] += vqx_[s][g][k];
                  }
                }
              }
            }

            // update digamma(vbeta)
            this.dg_vbtkc_[t][k][c] = MyGamma_digamma(vbeta_[t][k][c]);
            // dg_sum_vbtkc_[t][k] stores summand for every c
            this.dg_sum_vbtkc_[t][k] += vbeta_[t][k][c];
          }

          this.dg_sum_vbtkc_[t][k] = MyGamma_digamma(this.dg_sum_vbtkc_[t][k]);
        }
      }
      // M step by the Newton-Raphson method. (Maximize F(q) with respect to alpha and
      // beta.)
      // Newton-Raphson method for alpha
      if (ESTIMATE_ALPHA) {
        computeNewtonRaphson(alpha_, valpha_, dg_vask_, dg_sum_vask_);
      }

      // Newton-Raphson method for beta
      if (ESTIMATE_BETA) {
        for (int t = 0; t < T_; ++t) {
          computeNewtonRaphson(beta_[t], vbeta_[t], dg_vbtkc_[t], dg_sum_vbtkc_[t]);
        }
      }

      // compute the current value of F(q)
      double currentFq = computeFq();

      if (VBEM_MIN_ITERATION <= iteration) {
        if ((currentFq - prevFq) / Math.abs(prevFq) < THRESHOLD_FQ) {
          convergence = true;
        }

        if (iteration >= VBEM_MAX_ITERATION) {
          convergence = true;
        }
      }

      prevFq = currentFq;
      ++iteration;
    }

    return prevFq;
  }

  static private void computeNewtonRaphson(double[] a, final double[][] va, final double[][] dg_vask,
      final double[] dg_sum_vask) {
    int K = a.length;
    int S = va.length;
    double[] gradient = new double[K];

    // digamma of hyper-parameters alpha
    // dg_ak_[k] = digamma(alpha_k)
    double[] dg_ak = new double[K];
    // dg_sum_ak_ = digamma(sum_k alpha_k)
    double dg_sum_ak;

    // trigamma of hyper-parameters alpha and beta
    // dg_ak_[k] = trigamma(alpha_k)
    double[] tg_ak = new double[K];
    // dg_sum_ak_ = trigamma(sum_k alpha_k)
    double tg_sum_ak;

    int iteration = 0;
    double diffAlphaToCheckConvergence = 0;

    do {
      // precompute digamma and trigamma functions for hyper-parameters alpha and beta
      double sum_ak = 0;

      for (int k = 0; k < K; ++k) {
        sum_ak += a[k];
        dg_ak[k] = MyGamma_digamma(a[k]);
        tg_ak[k] = MyGamma_trigamma(a[k]);

        if (tg_ak[k] == 0) {
          System.err.println("trigamma=0, a[k]=" + a[k]);
        }
      }

      dg_sum_ak = MyGamma_digamma(sum_ak);
      tg_sum_ak = MyGamma_trigamma(sum_ak);

      // precompute gradient
      double gradient_common_term = S * dg_sum_ak;

      for (int s = 0; s < S; ++s) {
        gradient_common_term -= dg_sum_vask[s];
      }

      for (int k = 0; k < K; ++k) {
        gradient[k] = gradient_common_term - S * dg_ak[k];

        for (int s = 0; s < S; ++s) {
          gradient[k] += dg_vask[s][k];
        }
      }

      // compute common term (b)
      double c = 0;
      double d = 0;

      for (int k = 0; k < K; ++k) {
        d += 1 / (-S * tg_ak[k]);
        c += gradient[k] / (-S * tg_ak[k]);
      }

      double b = c / (1 / (S * tg_sum_ak) + d);

      // update alpha
      diffAlphaToCheckConvergence = 0;

      for (int k = 0; k < K; ++k) {
        double oldAlpha = a[k];
        // a[k] -= (gradient[k]-b)/(-S*tg_ak[k]);
        double direction = (gradient[k] - b) / (-S * tg_ak[k]);
        double newtonStepSize = 1;
        int i = 0;

        do {
          a[k] -= newtonStepSize * direction;
          newtonStepSize *= 0.5;
          ++i;
        } while (a[k] <= 0 && i <= 10);

        if (a[k] <= 0 || 10 <= a[k]) {
          a[k] = oldAlpha / 2.0;
        }

        diffAlphaToCheckConvergence += Math.abs(a[k] - oldAlpha) / (double) K;
      }

      ++iteration;
    } while (diffAlphaToCheckConvergence >= THRESHOLD_ALPHA && iteration < NEWTON_MAX_ITERATION);
  }

  private void printHyperParameters(PrintWriter pw) {
    pw.println("Alpha_k:");

    for (int k = 0; k < K_; ++k) {
      pw.println("\t" + alpha_[k]);
    }

    pw.println();
    pw.println("Beta_{t,c}:");

    for (int t = 0; t < T_; ++t) {
      pw.print("postWGD:" + (t + 1));

      for (int c = 0; c < Ct_[t]; ++c) {
        pw.print("\t" + beta_[t][c]);
      }

      pw.println();
    }
  }

  private void printMaxProbSegmentClusters(PrintWriter pw) {
    ArrayList<ArrayList<String>> segClusters = new ArrayList<ArrayList<String>>(K_);

    for (int i = 0; i < K_; ++i) {
      segClusters.add(new ArrayList<String>());
    }

    for (int s = 0; s < S_; ++s) {
      double max = -Double.MAX_VALUE;
      int argmaxk = -1;

      for (int k = 0; k < K_; ++k) {
        if (valpha_[s][k] >= max) {
          argmaxk = k;
          max = valpha_[s][k];
        }
      }

      Segment seg = this.segmentList_.get(s);
      segClusters.get(argmaxk).add(seg.getID());
    }
    for (int i = 0; i < K_; ++i) {
      ArrayList<String> cluster = segClusters.get(i);
      String sep = "";

      for (int j = 0; j < cluster.size(); ++j) {
        String segID = cluster.get(j);
        pw.print(sep + segID);
        sep = ",";
      }

      pw.println();
    }
  }

  private void printVariationalAlpha(PrintWriter pw) {
    for (int s = 0; s < S_; ++s) {
      Segment seg = this.segmentList_.get(s);
      pw.print(seg.getID());

      for (int k = 0; k < K_; ++k) {
        pw.print("\t" + valpha_[s][k]);
      }

      pw.println();
    }
  }

  private void printEU(PrintWriter pw) {
    for (int s = 0; s < S_; ++s) {
      Segment seg = this.segmentList_.get(s);
      pw.print(seg.getID());

      double sum = 0;

      for (int k = 0; k < K_; ++k) {
        sum += valpha_[s][k];
      }

      for (int k = 0; k < K_; ++k) {
        double p = valpha_[s][k] / sum;
        pw.print("\t" + p);
      }

      pw.println();
    }
  }

  private void printVariationalBeta(PrintWriter pw) {
    for (int t = 0; t < T_; ++t) {
      pw.println("postWGD species: " + (t + 1));

      for (int k = 0; k < K_; ++k) {
        pw.print("Anc" + (k + 1));

        for (int c = 0; c < Ct_[t]; ++c) {
          pw.print("\t" + vbeta_[t][k][c]);
        }

        pw.println();
      }

      pw.println();
    }
  }

  private void printV(PrintWriter pw) {
    for (int t = 0; t < T_; ++t) {
      pw.println("postWGD species: " + (t + 1));

      for (int k = 0; k < K_; ++k) {
        pw.print("Anc" + (k + 1));
        double sum = 0;

        for (int c = 0; c < Ct_[t]; ++c) {
          sum += vbeta_[t][k][c];
        }

        for (int c = 0; c < Ct_[t]; ++c) {
          double p = vbeta_[t][k][c] / sum;
          pw.print("\t" + p);
        }

        pw.println();
      }

      pw.println();
    }
  }

  private void printWordPosterior(PrintWriter pw) {
    for (int t = 0; t < T_; ++t) {
      pw.println("postWGD species " + (t + 1) + ":\t" + this.postWgdSpeciesNameList_[t]);
      double[][] prob = new double[K_][Ct_[t]];

      for (int k = 0; k < K_; ++k) {
        double total_Vbeta_c = 0;

        for (int c = 0; c < Ct_[t]; ++c) {
          total_Vbeta_c += vbeta_[t][k][c];
        }

        for (int c = 0; c < Ct_[t]; ++c) {
          prob[k][c] = vbeta_[t][k][c] / total_Vbeta_c;
        }
      }

      for (int c = 0; c < Ct_[t]; ++c) {
        double total_prob_k = 0;

        for (int k = 0; k < K_; ++k) {
          total_prob_k += prob[k][c];
        }

        pw.print(c);

        for (int k = 0; k < K_; ++k) {
          double posterior = prob[k][c] / total_prob_k;
          pw.print("\t" + posterior);
        }

        pw.println();
      }

      pw.println();
    }
  }

  private void printVariationalQX(PrintWriter pw) {
    for (int s = 0; s < S_; ++s) {
      Segment seg = this.segmentList_.get(s);
      pw.println(seg.getID());

      for (int g = 0; g < seg.getLength(); ++g) {
        Gene gene = seg.getGene(g);
        pw.print(gene.getID());

        for (int k = 0; k < K_; ++k) {
          if (LOGSUMEXP)
            pw.print("\t" + Math.exp(vqx_[s][g][k]));
          else
            pw.print("\t" + vqx_[s][g][k]);
        }

        pw.println();
      }

      pw.println();
    }
  }

  public void printResults(String outdir) throws IOException {
    String fileVA = outdir + File.separator + "est_varAsk.txt";
    String fileU = outdir + File.separator + "est_Usk.txt";
    String fileVB = outdir + File.separator + "est_varBtkc.txt";
    String fileV = outdir + File.separator + "est_Vtkc.txt";
    String fileX = outdir + File.separator + "est_varXsg.txt";
    String fileAB = outdir + File.separator + "est_AB.txt";
    String fileW = outdir + File.separator + "est_Wtk.txt";
    String fileMaxProbRec = outdir + File.separator + "rec_maxProb.txt";
    PrintWriter pwVA = new PrintWriter(new FileWriter(fileVA));
    PrintWriter pwU = new PrintWriter(new FileWriter(fileU));
    PrintWriter pwVB = new PrintWriter(new FileWriter(fileVB));
    PrintWriter pwV = new PrintWriter(new FileWriter(fileV));
    PrintWriter pwX = new PrintWriter(new FileWriter(fileX));
    PrintWriter pwAB = new PrintWriter(new FileWriter(fileAB));
    PrintWriter pwW = new PrintWriter(new FileWriter(fileW));
    PrintWriter pwMaxProbRec = new PrintWriter(new FileWriter(fileMaxProbRec));
    printVariationalAlpha(pwVA);
    printEU(pwU);
    printVariationalBeta(pwVB);
    printV(pwV);
    printWordPosterior(pwW);
    printVariationalQX(pwX);
    printHyperParameters(pwAB);
    printMaxProbSegmentClusters(pwMaxProbRec);
    pwVA.flush();
    pwVA.close();
    pwU.flush();
    pwU.close();
    pwVB.flush();
    pwVB.close();
    pwV.flush();
    pwV.close();
    pwX.flush();
    pwX.close();
    pwAB.flush();
    pwAB.close();
    pwW.flush();
    pwW.close();
    pwMaxProbRec.flush();
    pwMaxProbRec.close();
  }

  public int getGeneSize() {
    return this.numberOfGenes_;
  }

  public int getOrthologSize() {
    return this.numberOfOrthologs_;
  }

  static private void printUsage() {
    System.out.println("usage: ");
    System.out.println("  -?                   (=print usage)");
    System.out.println("  -i FILE              (=input synteny data file)");
    System.out.println("  -o DIR               (=output directory)");
    System.out.println("  -K INT               (=number of preWGD ancestral chromosomes)");
    System.out.println("  -T INT               (=number of postWGD species)");
    System.out.println("  -S STRING,...,STRING (=list of post-WGD species names t=1,...,T");
    System.out.println("  -C INT,...,INT       (=numbers of chromosomes in postWGD genomes t=1,...,T");
    System.out.println("  -D INT,...,INT       (=numbers of maximum co-orthologs in postWGD genomes t=1,...,T)");
    System.out.println("  -A DOUBLE            (=hyper-parameter for the symmetric dirichlet prior alpha)");
    System.out.println("  -B DOUBLE            (=hyper-parameter for the symmetric dirichlet prior beta)");
    System.out.println(
        "  -l DOUBLE            (=segments are excluded from shorter ones until the ratio of the total orthologs are excluded)");
    System.out
        .println("  -E BOOLEAN,BOOLEAN   (=Estimate hyper-parameter values by Newton-Raphson method. (Alpha,Beta))");
    System.out.println("  -I CLUSTERING/RANDOM_BETA/RANDOM_QX (=initialization method)");
    System.out.println("  -M VBEM/CVB/CVB0     (=inference method)");
  }

  public static void main(String[] args) throws IOException {
    ArrayList<String> inputSyntenyFileList = new ArrayList<String>();
    String outputDir = null;
    int K = 10;
    int T = 1;
    int[] Ct = null;
    int[] Dt = null;
    double alpha = 1;
    double beta = 1;
    double minSegLengthRatio = -Double.MAX_VALUE;
    String[] postWgdSpeciesNameList = null;
    boolean estimateAlpha = true;
    boolean estimateBeta = true;
    InitializationMethod initializationMethod = InitializationMethod.CLUSTERING;
    InferenceMethod inferenceMethod = InferenceMethod.CVB0;
    InitialSegmentClustering.DistanceMethod clusteringDistanceMethod = InitialSegmentClustering.DistanceMethod.SQUARE_ERROR;

    for (int i = 0; i < args.length; ++i) {
      switch (args[i].charAt(1)) {
        case 'i':
          inputSyntenyFileList.add(args[++i]);
          break;
        case 'o':
          outputDir = args[++i];
          break;
        case 'K':
          K = Integer.parseInt(args[++i]);
          break;
        case 'T':
          T = Integer.parseInt(args[++i]);
          break;
        case 'S':
          postWgdSpeciesNameList = args[++i].split(",");
          break;
        case 'C':
          String[] values = args[++i].split(",");
          Ct = new int[values.length];

          for (int j = 0; j < values.length; ++j) {
            Ct[j] = Integer.parseInt(values[j]);
          }

          break;
        case 'D':
          values = args[++i].split(",");
          Dt = new int[values.length];

          for (int j = 0; j < values.length; ++j) {
            Dt[j] = Integer.parseInt(values[j]);
          }

          break;
        case 'A':
          alpha = Double.parseDouble(args[++i]);
          break;
        case 'B':
          beta = Double.parseDouble(args[++i]);
          break;
        case 'E':
          values = args[++i].split(",");
          estimateAlpha = Boolean.valueOf(values[0]);
          estimateBeta = Boolean.valueOf(values[1]);
          break;
        case 'l':
          minSegLengthRatio = Double.parseDouble(args[++i]);
          break;
        case 'I':
          initializationMethod = InitializationMethod.valueOf(args[++i]);
          break;
        case 'M':
          inferenceMethod = InferenceMethod.valueOf(args[++i]);
          break;
        case 'X':
          clusteringDistanceMethod = InitialSegmentClustering.DistanceMethod.valueOf(args[++i]);
          break;
        case '?':
          printUsage();
          System.exit(0);
        default:
          System.err.println("Error: args[" + i + "]=" + args[i]);
          printUsage();
          System.exit(1);
      }
    }

    // print command line arguments
    if (VERBOSE) {
      System.out.print("Command line arguments:");

      for (int i = 0; i < args.length; ++i) {
        if (args[i].charAt(0) == '-') {
          System.out.println();
        } else {
          System.out.print(' ');
        }

        System.out.print(args[i]);
      }

      System.out.println();
    }

    if (inputSyntenyFileList.isEmpty() || outputDir == null || postWgdSpeciesNameList == null) {
      System.err.println("option error.");
      printUsage();

      System.exit(1);
    }

    File dir = new File(outputDir);

    MSynRec model = new MSynRec(postWgdSpeciesNameList, K, T, Ct, Dt);
    model.setWorkingDirectory(dir.getPath());
    model.setInitializationMethod(initializationMethod);
    model.setEstimateHyperParamters(estimateAlpha, estimateBeta);
    model.setInferenceMethod(inferenceMethod);

    for (String inputSyntenyFile : inputSyntenyFileList) {
      model.readSyntenyData(inputSyntenyFile);
    }

    if (minSegLengthRatio >= 0) {
      InitialSegmentClustering.setMinSegmentLengthRatio(minSegLengthRatio);
    }

    InitialSegmentClustering.setDistanceMethod(clusteringDistanceMethod);
    model.setPrior(alpha, beta);

    // double Fq = -1 + MIN_FQ_TO_AVOID_SUBOPTIMAL_INFERENCE;

    switch (inferenceMethod) {
      case CVB:
      case CVB0:
        model.computeCVB();
        break;
      case VBEM:
        model.computeVBEM();
        break;
    }

    model.printResults(dir.getPath());
  }

  static private double MyGamma_digamma(final double x) {
    if (Double.isNaN(x)) {
      throw new RuntimeException("digamma(NaN)");
    }

    if (x < -Double.MAX_VALUE) {
      System.err.println("digamma(" + x + ")");
      throw new RuntimeException();
    }

    double y = Gamma.digamma(x);

    return y;
  }

  static private double MyGamma_trigamma(final double x) {
    double y = Gamma.trigamma(x);

    return y;
  }

  static private double MyGamma_logGamma(final double x) {
    if (Double.isNaN(x)) {
      throw new RuntimeException("logGamma(NaN)");
    }

    double y = Gamma.logGamma(x);

    return y;
  }

  static private double logsumexp(double log_a, double log_b) {
    double logsumexp = 0;

    if (Double.isNaN(log_a) || Double.isInfinite(log_a))
      logsumexp = log_b;
    else if (Double.isNaN(log_b) || Double.isInfinite(log_b))
      logsumexp = log_a;
    else if (log_a > log_b) {
      logsumexp = log_a + Math.log(1 + Math.exp(log_b - log_a));
    } else {
      logsumexp = log_b + Math.log(1 + Math.exp(log_a - log_b));
    }

    return logsumexp;
  }
}
