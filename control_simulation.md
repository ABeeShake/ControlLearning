<script type="text/javascript"
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS_CHTML">
</script>
<script type="text/x-mathjax-config">
  MathJax.Hub.Config({
    tex2jax: {
      inlineMath: [['$','$'], ['\\(','\\)']],
      processEscapes: true},
      jax: ["input/TeX","input/MathML","input/AsciiMath","output/CommonHTML"],
      extensions: ["tex2jax.js","mml2jax.js","asciimath2jax.js","MathMenu.js","MathZoom.js","AssistiveMML.js", "[Contrib]/a11y/accessibility-menu.js"],
      TeX: {
      extensions: ["AMSmath.js","AMSsymbols.js","noErrors.js","noUndefined.js"],
      equationNumbers: {
      autoNumber: "AMS"
      }
    }
  });
</script>

$\renewcommand{\vec}[1]{\mathbf{#1}}$

# Control Simulation Notes

## Task 1:

Suppose we have a linear system

$$\dot{\vec{x}} = A\vec{x} + B\vec{u},$$

with states $\vec{x}$ and inputs $\vec{u}$. We will consider the case where $A$ and $B$ are unknown and must be estimated using data.

### Set Up:

1. Generate ground truth $A$ and $B$ with initial state $\vec{x}_0$ and pre-determined input sequence $\{\vec{u}\}_{t\in S}$, where $S$ is a set of points in time.
2. For a linear system, we use the controller $\vec{u}_t = -K\vec{x}_t$
3. Discretize the linear system as follows:

$$
\begin{align*}
\frac{\vec{x}_{t+\Delta t}- \vec{x}_{t}}{\Delta t} &= A\vec{x}_t + B\vec{u}_t\\
\implies \vec{x}_{t+\Delta t}- \vec{x}_{t} &= A\vec{x}_t\Delta t + B\vec{u}_t\Delta t\\
\implies \vec{x}_{t+\Delta t}&= (\Delta tA+I)\vec{x}_t + \Delta tB\vec{u}_t\\
\implies \vec{x}_{t+\Delta t}&= (\Delta tA+I -\Delta tBK)\vec{x}_t\\
\end{align*}
$$
4. Use the elements from (1.) to generate sequences $\{\vec{x}_t\}_{t\in S}$.
5. *Optional*: Add noise when generating $\{\vec{x}_t,\vec{u}_t\}_{t\in S}$.

### Estimation Process:

