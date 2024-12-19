## Overview

R code used for the machine learning portion of: 

Kieran TJ, Sun X, Maines TR, Belser JA. Optimal thresholds and key parameters for predicting influenza A virus transmission events in ferrets. npj Viruses 2, 64 (2024). https://doi.org/10.1038/s44298-024-00074-w

This project also includes CSV files for summary data & metrics used to create figures from the manuscript in the "ResultsSummary_Rinputs" directory. 

Study makes use of previously published data:
## References and Resources
Kieran TJ, Sun X, Maines TR, Belser JA. Machine learning approaches for influenza A virus risk assessment identifies predictive correlates using ferret model in vivo data. Commun Biol 7, 927 (2024). https://doi.org/10.1038/s42003-024-06629-0

https://github.com/CDCgov/machine-learning-influenza-ferret-model/tree/main

Kieran TJ, Sun X, Creager HM, Tumpey TM, Maines TR, Belser JA. An aggregated dataset of serial morbidity and titer measurements from influenza A virus-infected ferrets. Sci Data 11, 510 (2024). https://doi.org/10.1038/s41597-024-03256-6

Kieran TJ, Sun X, Maines TR, Beauchemin CAA, Belser JA. Exploring associations between viral titer measurements and disease outcomes in ferrets inoculated with 125 contemporary influenza A viruses. J Virol 98:e01661-23. (2024). (PMID 38240592)
([https://doi.org/10.1128/jvi.01661-23](https://doi.org/10.1128/jvi.01661-23))

## Note
R script uses predicted sialic acid binding preference (RBS) and predicted polymerase activity (PBS) markers (based on molecular sequence), along with selected HA and PB2 gene markers, that are not included in the CSV file from Sci Data, but may be cross-referenced as reported in Kieran et al (PMID 38240592).

## Manuscript Abstract
Although assessments of influenza A virus transmissibility in the ferret model play a critical role in pandemic risk evaluations, few studies have investigated which virological data collected from virus-inoculated animals are most predictive of subsequent virus transmission to naïve contacts. We compiled viral titer data from >475 ferrets inoculated with 97 contemporary IAV (including high- and low-pathogenicity avian, swine-origin, and human viruses of multiple HA subtypes) that served as donors for assessments of virus transmission in the presence of direct contact (DCT) or via respiratory droplets (RDT). A diversity of molecular determinants, clinical parameters, and infectious titer measurements and derived quantities were examined to identify which metrics were most statistically supported with transmission outcome. Higher viral loads in nasal wash (NW) specimens were strongly associated with higher transmission frequencies in DCT, but not RDT models. However, viruses that reached peak titers in NW specimens early (day 1 p.i.) were strongly associated with higher transmission in both models. Interestingly, viruses with ‘intermediate’ transmission outcomes (33-66%) had NW titers and derived quantities more similar to non-transmissible viruses (<33%) in a DCT setting, but with efficiently transmissible viruses (>67%) in a RDT setting. Machine learning was employed to further assess the predictive role of summary measures and varied interpretation of intermediate transmission outcomes in both DCT and RDT models, with models employing these different thresholds yielding high performance metrics against both internal and external datasets. Collectively, these findings suggest that higher viral load in inoculated animals can be predictive of DCT outcomes, whereas the timing of when peak titers are detected in inoculated animals can inform RDT outcomes. Identification that intermediate transmission outcomes should be contextualized relative to the transmission mode assessed provides needed refinement towards improving interpretation of ferret transmission studies in the context of pandemic risk assessment.

##
##
##
  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).
