package math;

public class SpecialFunctions {
  private static final double HALF_LOG_2_PI = 0.5 * Math.log(2.0 * Math.PI);

  public static final double LANCZOS_G = 607.0 / 128.0;

  private static final double[] LANCZOS = {
      0.99999999999999709182,
      57.156235665862923517,
      -59.597960355475491248,
      14.136097974741747174,
      -0.49191381609762019978,
      .33994649984811888699e-4,
      .46523628927048575665e-4,
      -.98374475304879564677e-4,
      .15808870322491248884e-3,
      -.21026444172410488319e-3,
      .21743961811521264320e-3,
      -.16431810653676389022e-3,
      .84418223983852743293e-4,
      -.26190838401581408670e-4,
      .36899182659531622704e-5,
  };

  // reference: http://home.apache.org/~luc/commons-math-3.6-RC2-site/jacoco/org.apache.commons.math3.special/Gamma.java.html
  public static double lanczos(final double x) {
    double sum = 0;

    for (int i = LANCZOS.length - 1; i > 0; --i) {
      sum += LANCZOS[i] / (x + i);
    }

    return sum + LANCZOS[0];
  }

  // NOTE: x will always be large
  // reference: http://home.apache.org/~luc/commons-math-3.6-RC2-site/jacoco/org.apache.commons.math3.special/Gamma.java.html
  public static double logGamma(double x) {
    double sum = lanczos(x);
    double tmp = x + LANCZOS_G + 0.5;

    return ((x + 0.5) * Math.log(tmp)) - tmp + HALF_LOG_2_PI + Math.log(sum / x);
  }

  // reference: https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
  public static double logsumexp(double a, double b) {
    double c = a > b ? a : b;
    double sum = Math.exp(a - c) + Math.exp(b - c);

    return Math.log(sum) + c;
  }

  public static double logFactorial(int d) {
    return logGamma(d + 1);
  }

  public static double logCombination(int N, int n) {
    if (n < 0 || N < n) {
      return Double.NEGATIVE_INFINITY;
    }

    return logFactorial(N) - logFactorial(n) - logFactorial(N - n);
  }

  public static double logHypergeometricUpperProbability(int n, int m, int N, int i) {
    double logsump = Double.NEGATIVE_INFINITY;
    double constant = logCombination(n + m, N);

    for (int j = i; j <= n; ++j) {
      double logp = 0;
      logp += logCombination(n, j);
      logp += logCombination(m, N - j);
      logp -= constant;

      logsump = logsumexp(logsump, logp);
    }

    return logsump;
  }

  public static double logHypergeometricLowerProbability(int n, int m, int N, int i) {
    double logsump = Double.NEGATIVE_INFINITY;
    double constant = logCombination(n + m, N);

    for (int j = 0; j <= i; ++j) {
      double logp = 0;
      logp += logCombination(n, j);
      logp += logCombination(m, N - j);
      logp -= constant;

      logsump = logsumexp(logsump, logp);
    }

    return logsump;
  }
}
