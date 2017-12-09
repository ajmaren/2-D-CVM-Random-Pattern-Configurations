# 2-D-CVM-Random-Pattern-Configurations
Code, images, PPTs, spreadsheets, and PDF documents to create randomly-generated toy problem arrays (16x16, or 256 units) for 2-D CVM analysis. CVM stands for "Cluster Variation Method," initially developed by R. Kikuchi. 
The 2-D CVM is a two-dimensional array of on/off (A/B, or 1/0) units. They are arranged in a zigzag grid. To see examples, please consult the documents listed at the end of this ReadMe. 
It is possible to compute thermodynamic variables, such as entropy, enthalpy, and free energy for a CVM system. The key difference between the CVM and "normal" Ising models is that in the CVM, the entropy is more complex. Instead of just being a function of the proportional number of units in states A and B (x1 and x2), it is also a function of: 
- nearest neighbors, or A-A, A-B / B-A, and B-B pairs
- next-nearest neighbors, or A--A, A--B / B--A, and B--B pairs, and
- triplets, such as A-B-A, B-A-A, etc. (There are six kinds of triplet patterns.) 
Not surprisingly, the entropy term is much more complex than usual. 

I have developed an analytic solution for the entropy at the equiprobable distribution points (x1 = x2 = 0.5). Unfortunately, there are difficulties with systems that have equiprobable distributions. (Details forthcoming in various papers.) 

For all other distributions, it is more complex to determine the configuration variables (the nearest neighbors, next-nearest neighbors, and triplets) that correspond to specific enthalpy parameters and to free energy minima. What is important is that, if the system is not at an equiprobable distribution (x1 = x2 = 0.5), the point at which a system is at free energy minimum must be found computationally. Further, this is not as straightforward as solving for a minimum in a simple equation; e.g., a Newton-Raphson method would not work.

In fact, finding good solutions require dealing with network topographies. 

This GitHub repository contains useful code and documentation. 

If you use this material, please credit me: Alianna J. Maren.

For the latest in results, please go to: --- www.aliannajmaren.com --- and Opt-In using the Opt-In form in the right-hand sidebar. 
You can also see various results posted in the blog; they can be found by using the blog categories. 

For more information or to discuss the CVM methods, reach me at:
-- alianna@aliannajmaren.com
-- alianna.maren@northwestern.edu


Online Journal Article: 
- Maren, A.J. The Cluster Variation Method: A Primer for Neuroscientists, Brain Sci. 2016, 6(4), 44; doi:10.3390/brainsci6040044, http://www.mdpi.com/2076-3425/6/4/44
