# Strobl et al (2020). Turnover modulates the need for a cost of resistance in adaptive therapy

<p align="center">
  <img src="https://cancerres.aacrjournals.org/sites/default/files/styles/medium/public/highwire/canres/81/4.cover-source.jpg", alt="Cover of the print-version of the edition of Cancer Research, featuring a collage of phase planes from our paper."/>
</p>

This repository contains the code to generate the figures in our research paper: Strobl et al (2020), [*Turnover modulates the need for a cost of resistance in adaptive therapy*](https://cancerres.aacrjournals.org/content/81/4/1135) published in *Cancer Research*. A pre-print of this paper is available [here](https://www.biorxiv.org/content/biorxiv/early/2020/01/29/2020.01.22.914366.full.pdf). 

For each figure, as well as for the ingredients for the cover image, there is a dedicated Jupyter notebook. The files in the `utils` folder contain all the main functions for simulating the model and plotting its results. For reproducibility, we have also posted here the version of the clinical data by Bruchovsky et al (2007) which we downloaded from [here](http://www.nicholasbruchovsky.com/clinicalResearch.html) and analysed in July 2020. This is located in `data/clinicalData`. We thank the authors for making this data publicly available.

To run run the code, install the required dependencies via
`pip install -r requirements.txt`. Then start a jupyter notebook and fire away.

#### References

- ﻿Bruchovsky, N., Klotz, L., Crook, J., Malone, S., Ludgate, C., Morris, W. J., … Goldenberg, S. L. (2006). Final results of the Canadian prospective Phase II trial of intermittent androgen suppression for men in biochemical recurrence after radiotherapy for locally advanced prostate cancer: Clinical parameters. Cancer, 107(2), 389–395. https://doi.org/10.1002/cncr.21989
- ﻿Strobl, M. A. R., West, J., Viossat, Y., Damaghi, M., Robertson-Tessi, M., Brown, J. S., Gatenby, R. A., Maini, P. K., Anderson, A. R. A. (2020). Turnover modulates the need for a cost of resistance in adaptive therapy. Cancer Research, (OnlineFirst), November 10, 2020. https://doi.org/10.1158/0008-5472.can-20-0806