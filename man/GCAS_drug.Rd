<<<<<<< HEAD
\name{GCAS_drug}
\alias{GCAS_drug}
\docType{data}
\title{
Anti-tumor drug sensitivity of all samples
}
\description{
Anti-tumor drug sensitivity of all samples calculated based on oncoPredict package and GDSC database.
}
\usage{data("GCAS_drug")}
\format{
  A data frame with 19816 observations on the following 199 variables.
  \describe{
    \item{\code{ID}}{a character vector}
    \item{\code{Camptothecin_1003}}{a numeric vector}
    \item{\code{Vinblastine_1004}}{a numeric vector}
    \item{\code{Cisplatin_1005}}{a numeric vector}
    \item{\code{Cytarabine_1006}}{a numeric vector}
    \item{\code{Docetaxel_1007}}{a numeric vector}
    \item{\code{Gefitinib_1010}}{a numeric vector}
    \item{\code{Navitoclax_1011}}{a numeric vector}
    \item{\code{Vorinostat_1012}}{a numeric vector}
    \item{\code{Nilotinib_1013}}{a numeric vector}
    \item{\code{Olaparib_1017}}{a numeric vector}
    \item{\code{Axitinib_1021}}{a numeric vector}
    \item{\code{AZD7762_1022}}{a numeric vector}
    \item{\code{SB216763_1025}}{a numeric vector}
    \item{\samp{KU-55933_1030}}{a numeric vector}
    \item{\code{Afatinib_1032}}{a numeric vector}
    \item{\code{Staurosporine_1034}}{a numeric vector}
    \item{\samp{PLX-4720_1036}}{a numeric vector}
    \item{\code{NU7441_1038}}{a numeric vector}
    \item{\code{Doramapimod_1042}}{a numeric vector}
    \item{\samp{Wee1 Inhibitor_1046}}{a numeric vector}
    \item{\samp{Nutlin-3a (-)_1047}}{a numeric vector}
    \item{\code{Mirin_1048}}{a numeric vector}
    \item{\code{PD173074_1049}}{a numeric vector}
    \item{\code{ZM447439_1050}}{a numeric vector}
    \item{\code{Alisertib_1051}}{a numeric vector}
    \item{\samp{RO-3306_1052}}{a numeric vector}
    \item{\samp{MK-2206_1053}}{a numeric vector}
    \item{\code{Palbociclib_1054}}{a numeric vector}
    \item{\code{Dactolisib_1057}}{a numeric vector}
    \item{\code{Pictilisib_1058}}{a numeric vector}
    \item{\code{AZD8055_1059}}{a numeric vector}
    \item{\code{PD0325901_1060}}{a numeric vector}
    \item{\samp{Obatoclax Mesylate_1068}}{a numeric vector}
    \item{\samp{5-Fluorouracil_1073}}{a numeric vector}
    \item{\code{Dasatinib_1079}}{a numeric vector}
    \item{\code{Paclitaxel_1080}}{a numeric vector}
    \item{\code{Crizotinib_1083}}{a numeric vector}
    \item{\code{Rapamycin_1084}}{a numeric vector}
    \item{\code{Sorafenib_1085}}{a numeric vector}
    \item{\samp{BI-2536_1086}}{a numeric vector}
    \item{\code{Irinotecan_1088}}{a numeric vector}
    \item{\code{Oxaliplatin_1089}}{a numeric vector}
    \item{\samp{BMS-536924_1091}}{a numeric vector}
    \item{\code{GSK1904529A_1093}}{a numeric vector}
    \item{\code{Tozasertib_1096}}{a numeric vector}
    \item{\samp{PF-4708671_1129}}{a numeric vector}
    \item{\samp{PRIMA-1MET_1131}}{a numeric vector}
    \item{\code{Erlotinib_1168}}{a numeric vector}
    \item{\code{Niraparib_1177}}{a numeric vector}
    \item{\samp{MK-1775_1179}}{a numeric vector}
    \item{\code{Dinaciclib_1180}}{a numeric vector}
    \item{\code{Gemcitabine_1190}}{a numeric vector}
    \item{\code{Bortezomib_1191}}{a numeric vector}
    \item{\code{GSK269962A_1192}}{a numeric vector}
    \item{\code{SB505124_1194}}{a numeric vector}
    \item{\code{Tamoxifen_1199}}{a numeric vector}
    \item{\code{Fulvestrant_1200}}{a numeric vector}
    \item{\code{EPZ004777_1237}}{a numeric vector}
    \item{\samp{YK-4-279_1239}}{a numeric vector}
    \item{\code{Daporinad_1248}}{a numeric vector}
    \item{\samp{BMS-345541_1249}}{a numeric vector}
    \item{\code{AZ960_1250}}{a numeric vector}
    \item{\code{Talazoparib_1259}}{a numeric vector}
    \item{\code{XAV939_1268}}{a numeric vector}
    \item{\code{Trametinib_1372}}{a numeric vector}
    \item{\code{Dabrafenib_1373}}{a numeric vector}
    \item{\code{Temozolomide_1375}}{a numeric vector}
    \item{\code{AZD5438_1401}}{a numeric vector}
    \item{\code{IAP_5620_1428}}{a numeric vector}
    \item{\code{AZD2014_1441}}{a numeric vector}
    \item{\code{AZD1208_1449}}{a numeric vector}
    \item{\code{AZD1332_1463}}{a numeric vector}
    \item{\code{Ruxolitinib_1507}}{a numeric vector}
    \item{\code{Linsitinib_1510}}{a numeric vector}
    \item{\code{Epirubicin_1511}}{a numeric vector}
    \item{\code{Cyclophosphamide_1512}}{a numeric vector}
    \item{\code{Pevonedistat_1529}}{a numeric vector}
    \item{\code{Sapitinib_1549}}{a numeric vector}
    \item{\code{Uprosertib_1553}}{a numeric vector}
    \item{\code{LCL161_1557}}{a numeric vector}
    \item{\code{Lapatinib_1558}}{a numeric vector}
    \item{\code{Luminespib_1559}}{a numeric vector}
    \item{\code{Alpelisib_1560}}{a numeric vector}
    \item{\code{Taselisib_1561}}{a numeric vector}
    \item{\code{EPZ5676_1563}}{a numeric vector}
    \item{\code{SCH772984_1564}}{a numeric vector}
    \item{\samp{IWP-2_1576}}{a numeric vector}
    \item{\code{Leflunomide_1578}}{a numeric vector}
    \item{\code{Entinostat_1593}}{a numeric vector}
    \item{\samp{OSI-027_1594}}{a numeric vector}
    \item{\code{LGK974_1598}}{a numeric vector}
    \item{\samp{VE-822_1613}}{a numeric vector}
    \item{\code{WZ4003_1614}}{a numeric vector}
    \item{\code{CZC24832_1615}}{a numeric vector}
    \item{\code{AZD5582_1617}}{a numeric vector}
    \item{\code{GSK2606414_1618}}{a numeric vector}
    \item{\code{PFI3_1620}}{a numeric vector}
    \item{\samp{PCI-34051_1621}}{a numeric vector}
    \item{\samp{Wnt-C59_1622}}{a numeric vector}
    \item{\samp{I-BET-762_1624}}{a numeric vector}
    \item{\samp{RVX-208_1625}}{a numeric vector}
    \item{\code{OTX015_1626}}{a numeric vector}
    \item{\code{GSK343_1627}}{a numeric vector}
    \item{\code{ML323_1629}}{a numeric vector}
    \item{\code{Entospletinib_1630}}{a numeric vector}
    \item{\code{PRT062607_1631}}{a numeric vector}
    \item{\code{Ribociclib_1632}}{a numeric vector}
    \item{\samp{AGI-6780_1634}}{a numeric vector}
    \item{\samp{Picolinici-acid_1635}}{a numeric vector}
    \item{\code{AZD5153_1706}}{a numeric vector}
    \item{\code{CDK9_5576_1708}}{a numeric vector}
    \item{\code{CDK9_5038_1709}}{a numeric vector}
    \item{\code{Eg5_9814_1712}}{a numeric vector}
    \item{\code{ERK_2440_1713}}{a numeric vector}
    \item{\code{ERK_6604_1714}}{a numeric vector}
    \item{\code{IRAK4_4710_1716}}{a numeric vector}
    \item{\code{JAK1_8709_1718}}{a numeric vector}
    \item{\code{AZD5991_1720}}{a numeric vector}
    \item{\code{PAK_5339_1730}}{a numeric vector}
    \item{\code{TAF1_5496_1732}}{a numeric vector}
    \item{\code{ULK1_4989_1733}}{a numeric vector}
    \item{\code{VSP34_8731_1734}}{a numeric vector}
    \item{\code{Selumetinib_1736}}{a numeric vector}
    \item{\code{IGF1R_3801_1738}}{a numeric vector}
    \item{\code{JAK_8517_1739}}{a numeric vector}
    \item{\code{AZD4547_1786}}{a numeric vector}
    \item{\code{Ibrutinib_1799}}{a numeric vector}
    \item{\code{Zoledronate_1802}}{a numeric vector}
    \item{\code{Acetalax_1804}}{a numeric vector}
    \item{\code{Oxaliplatin_1806}}{a numeric vector}
    \item{\code{Carmustine_1807}}{a numeric vector}
    \item{\code{Topotecan_1808}}{a numeric vector}
    \item{\code{Teniposide_1809}}{a numeric vector}
    \item{\code{Mitoxantrone_1810}}{a numeric vector}
    \item{\code{Dactinomycin_1811}}{a numeric vector}
    \item{\code{Fludarabine_1813}}{a numeric vector}
    \item{\code{Nelarabine_1814}}{a numeric vector}
    \item{\code{Fulvestrant_1816}}{a numeric vector}
    \item{\code{Vincristine_1818}}{a numeric vector}
    \item{\code{Docetaxel_1819}}{a numeric vector}
    \item{\samp{Podophyllotoxin bromide_1825}}{a numeric vector}
    \item{\code{Dihydrorotenone_1827}}{a numeric vector}
    \item{\code{Gallibiscoquinazole_1830}}{a numeric vector}
    \item{\code{Elephantin_1835}}{a numeric vector}
    \item{\code{Sinularin_1838}}{a numeric vector}
    \item{\code{Sabutoclax_1849}}{a numeric vector}
    \item{\code{LY2109761_1852}}{a numeric vector}
    \item{\samp{OF-1_1853}}{a numeric vector}
    \item{\samp{MN-64_1854}}{a numeric vector}
    \item{\samp{KRAS (G12C) Inhibitor-12_1855}}{a numeric vector}
    \item{\samp{MG-132_1862}}{a numeric vector}
    \item{\samp{BDP-00009066_1866}}{a numeric vector}
    \item{\code{Buparlisib_1873}}{a numeric vector}
    \item{\code{Ulixertinib_1908}}{a numeric vector}
    \item{\code{Venetoclax_1909}}{a numeric vector}
    \item{\code{ABT737_1910}}{a numeric vector}
    \item{\code{Dactinomycin_1911}}{a numeric vector}
    \item{\code{Afuresertib_1912}}{a numeric vector}
    \item{\samp{AGI-5198_1913}}{a numeric vector}
    \item{\code{AZD3759_1915}}{a numeric vector}
    \item{\code{AZD5363_1916}}{a numeric vector}
    \item{\code{AZD6738_1917}}{a numeric vector}
    \item{\code{AZD8186_1918}}{a numeric vector}
    \item{\code{Osimertinib_1919}}{a numeric vector}
    \item{\code{Cediranib_1922}}{a numeric vector}
    \item{\code{Ipatasertib_1924}}{a numeric vector}
    \item{\code{GDC0810_1925}}{a numeric vector}
    \item{\samp{GNE-317_1926}}{a numeric vector}
    \item{\code{GSK2578215A_1927}}{a numeric vector}
    \item{\samp{I-BRD9_1928}}{a numeric vector}
    \item{\samp{Telomerase Inhibitor IX_1930}}{a numeric vector}
    \item{\samp{MIRA-1_1931}}{a numeric vector}
    \item{\samp{NVP-ADW742_1932}}{a numeric vector}
    \item{\code{P22077_1933}}{a numeric vector}
    \item{\code{Savolitinib_1936}}{a numeric vector}
    \item{\samp{UMI-77_1939}}{a numeric vector}
    \item{\code{WIKI4_1940}}{a numeric vector}
    \item{\samp{Sepantronium bromide_1941}}{a numeric vector}
    \item{\code{MIM1_1996}}{a numeric vector}
    \item{\samp{WEHI-539_1997}}{a numeric vector}
    \item{\samp{BPD-00008900_1998}}{a numeric vector}
    \item{\code{Foretinib_2040}}{a numeric vector}
    \item{\samp{BIBR-1532_2043}}{a numeric vector}
    \item{\code{Pyridostatin_2044}}{a numeric vector}
    \item{\samp{AMG-319_2045}}{a numeric vector}
    \item{\samp{MK-8776_2046}}{a numeric vector}
    \item{\code{Ulixertinib_2047}}{a numeric vector}
    \item{\code{Vinorelbine_2048}}{a numeric vector}
    \item{\samp{VX-11e_2096}}{a numeric vector}
    \item{\code{Uprosertib_2106}}{a numeric vector}
    \item{\code{LJI308_2107}}{a numeric vector}
    \item{\code{AZ6102_2109}}{a numeric vector}
    \item{\code{GSK591_2110}}{a numeric vector}
    \item{\code{VE821_2111}}{a numeric vector}
    \item{\code{AZD6482_2169}}{a numeric vector}
    \item{\code{AT13148_2170}}{a numeric vector}
    \item{\samp{BMS-754807_2171}}{a numeric vector}
    \item{\code{JQ1_2172}}{a numeric vector}
  }
}
\references{
Maeser, D., Gruener, R. F., & Huang, R. S. (2021). oncoPredict: an R package for predicting in vivo or cancer patient drug response and biomarkers from cell line screening data. Briefings in bioinformatics, 22(6), bbab260. https://doi.org/10.1093/bib/bbab260
}
\examples{
data(GCAS_drug)
str(GCAS_drug)
}
=======
\name{GCAS_drug}
\alias{GCAS_drug}
\docType{data}
\title{
Anti-tumor drug sensitivity of all samples
}
\description{
Anti-tumor drug sensitivity of all samples calculated based on oncoPredict package and GDSC database.
}
\usage{data("GCAS_drug")}
\format{
  A data frame with 19816 observations on the following 199 variables.
  \describe{
    \item{\code{ID}}{a character vector}
    \item{\code{Camptothecin_1003}}{a numeric vector}
    \item{\code{Vinblastine_1004}}{a numeric vector}
    \item{\code{Cisplatin_1005}}{a numeric vector}
    \item{\code{Cytarabine_1006}}{a numeric vector}
    \item{\code{Docetaxel_1007}}{a numeric vector}
    \item{\code{Gefitinib_1010}}{a numeric vector}
    \item{\code{Navitoclax_1011}}{a numeric vector}
    \item{\code{Vorinostat_1012}}{a numeric vector}
    \item{\code{Nilotinib_1013}}{a numeric vector}
    \item{\code{Olaparib_1017}}{a numeric vector}
    \item{\code{Axitinib_1021}}{a numeric vector}
    \item{\code{AZD7762_1022}}{a numeric vector}
    \item{\code{SB216763_1025}}{a numeric vector}
    \item{\samp{KU-55933_1030}}{a numeric vector}
    \item{\code{Afatinib_1032}}{a numeric vector}
    \item{\code{Staurosporine_1034}}{a numeric vector}
    \item{\samp{PLX-4720_1036}}{a numeric vector}
    \item{\code{NU7441_1038}}{a numeric vector}
    \item{\code{Doramapimod_1042}}{a numeric vector}
    \item{\samp{Wee1 Inhibitor_1046}}{a numeric vector}
    \item{\samp{Nutlin-3a (-)_1047}}{a numeric vector}
    \item{\code{Mirin_1048}}{a numeric vector}
    \item{\code{PD173074_1049}}{a numeric vector}
    \item{\code{ZM447439_1050}}{a numeric vector}
    \item{\code{Alisertib_1051}}{a numeric vector}
    \item{\samp{RO-3306_1052}}{a numeric vector}
    \item{\samp{MK-2206_1053}}{a numeric vector}
    \item{\code{Palbociclib_1054}}{a numeric vector}
    \item{\code{Dactolisib_1057}}{a numeric vector}
    \item{\code{Pictilisib_1058}}{a numeric vector}
    \item{\code{AZD8055_1059}}{a numeric vector}
    \item{\code{PD0325901_1060}}{a numeric vector}
    \item{\samp{Obatoclax Mesylate_1068}}{a numeric vector}
    \item{\samp{5-Fluorouracil_1073}}{a numeric vector}
    \item{\code{Dasatinib_1079}}{a numeric vector}
    \item{\code{Paclitaxel_1080}}{a numeric vector}
    \item{\code{Crizotinib_1083}}{a numeric vector}
    \item{\code{Rapamycin_1084}}{a numeric vector}
    \item{\code{Sorafenib_1085}}{a numeric vector}
    \item{\samp{BI-2536_1086}}{a numeric vector}
    \item{\code{Irinotecan_1088}}{a numeric vector}
    \item{\code{Oxaliplatin_1089}}{a numeric vector}
    \item{\samp{BMS-536924_1091}}{a numeric vector}
    \item{\code{GSK1904529A_1093}}{a numeric vector}
    \item{\code{Tozasertib_1096}}{a numeric vector}
    \item{\samp{PF-4708671_1129}}{a numeric vector}
    \item{\samp{PRIMA-1MET_1131}}{a numeric vector}
    \item{\code{Erlotinib_1168}}{a numeric vector}
    \item{\code{Niraparib_1177}}{a numeric vector}
    \item{\samp{MK-1775_1179}}{a numeric vector}
    \item{\code{Dinaciclib_1180}}{a numeric vector}
    \item{\code{Gemcitabine_1190}}{a numeric vector}
    \item{\code{Bortezomib_1191}}{a numeric vector}
    \item{\code{GSK269962A_1192}}{a numeric vector}
    \item{\code{SB505124_1194}}{a numeric vector}
    \item{\code{Tamoxifen_1199}}{a numeric vector}
    \item{\code{Fulvestrant_1200}}{a numeric vector}
    \item{\code{EPZ004777_1237}}{a numeric vector}
    \item{\samp{YK-4-279_1239}}{a numeric vector}
    \item{\code{Daporinad_1248}}{a numeric vector}
    \item{\samp{BMS-345541_1249}}{a numeric vector}
    \item{\code{AZ960_1250}}{a numeric vector}
    \item{\code{Talazoparib_1259}}{a numeric vector}
    \item{\code{XAV939_1268}}{a numeric vector}
    \item{\code{Trametinib_1372}}{a numeric vector}
    \item{\code{Dabrafenib_1373}}{a numeric vector}
    \item{\code{Temozolomide_1375}}{a numeric vector}
    \item{\code{AZD5438_1401}}{a numeric vector}
    \item{\code{IAP_5620_1428}}{a numeric vector}
    \item{\code{AZD2014_1441}}{a numeric vector}
    \item{\code{AZD1208_1449}}{a numeric vector}
    \item{\code{AZD1332_1463}}{a numeric vector}
    \item{\code{Ruxolitinib_1507}}{a numeric vector}
    \item{\code{Linsitinib_1510}}{a numeric vector}
    \item{\code{Epirubicin_1511}}{a numeric vector}
    \item{\code{Cyclophosphamide_1512}}{a numeric vector}
    \item{\code{Pevonedistat_1529}}{a numeric vector}
    \item{\code{Sapitinib_1549}}{a numeric vector}
    \item{\code{Uprosertib_1553}}{a numeric vector}
    \item{\code{LCL161_1557}}{a numeric vector}
    \item{\code{Lapatinib_1558}}{a numeric vector}
    \item{\code{Luminespib_1559}}{a numeric vector}
    \item{\code{Alpelisib_1560}}{a numeric vector}
    \item{\code{Taselisib_1561}}{a numeric vector}
    \item{\code{EPZ5676_1563}}{a numeric vector}
    \item{\code{SCH772984_1564}}{a numeric vector}
    \item{\samp{IWP-2_1576}}{a numeric vector}
    \item{\code{Leflunomide_1578}}{a numeric vector}
    \item{\code{Entinostat_1593}}{a numeric vector}
    \item{\samp{OSI-027_1594}}{a numeric vector}
    \item{\code{LGK974_1598}}{a numeric vector}
    \item{\samp{VE-822_1613}}{a numeric vector}
    \item{\code{WZ4003_1614}}{a numeric vector}
    \item{\code{CZC24832_1615}}{a numeric vector}
    \item{\code{AZD5582_1617}}{a numeric vector}
    \item{\code{GSK2606414_1618}}{a numeric vector}
    \item{\code{PFI3_1620}}{a numeric vector}
    \item{\samp{PCI-34051_1621}}{a numeric vector}
    \item{\samp{Wnt-C59_1622}}{a numeric vector}
    \item{\samp{I-BET-762_1624}}{a numeric vector}
    \item{\samp{RVX-208_1625}}{a numeric vector}
    \item{\code{OTX015_1626}}{a numeric vector}
    \item{\code{GSK343_1627}}{a numeric vector}
    \item{\code{ML323_1629}}{a numeric vector}
    \item{\code{Entospletinib_1630}}{a numeric vector}
    \item{\code{PRT062607_1631}}{a numeric vector}
    \item{\code{Ribociclib_1632}}{a numeric vector}
    \item{\samp{AGI-6780_1634}}{a numeric vector}
    \item{\samp{Picolinici-acid_1635}}{a numeric vector}
    \item{\code{AZD5153_1706}}{a numeric vector}
    \item{\code{CDK9_5576_1708}}{a numeric vector}
    \item{\code{CDK9_5038_1709}}{a numeric vector}
    \item{\code{Eg5_9814_1712}}{a numeric vector}
    \item{\code{ERK_2440_1713}}{a numeric vector}
    \item{\code{ERK_6604_1714}}{a numeric vector}
    \item{\code{IRAK4_4710_1716}}{a numeric vector}
    \item{\code{JAK1_8709_1718}}{a numeric vector}
    \item{\code{AZD5991_1720}}{a numeric vector}
    \item{\code{PAK_5339_1730}}{a numeric vector}
    \item{\code{TAF1_5496_1732}}{a numeric vector}
    \item{\code{ULK1_4989_1733}}{a numeric vector}
    \item{\code{VSP34_8731_1734}}{a numeric vector}
    \item{\code{Selumetinib_1736}}{a numeric vector}
    \item{\code{IGF1R_3801_1738}}{a numeric vector}
    \item{\code{JAK_8517_1739}}{a numeric vector}
    \item{\code{AZD4547_1786}}{a numeric vector}
    \item{\code{Ibrutinib_1799}}{a numeric vector}
    \item{\code{Zoledronate_1802}}{a numeric vector}
    \item{\code{Acetalax_1804}}{a numeric vector}
    \item{\code{Oxaliplatin_1806}}{a numeric vector}
    \item{\code{Carmustine_1807}}{a numeric vector}
    \item{\code{Topotecan_1808}}{a numeric vector}
    \item{\code{Teniposide_1809}}{a numeric vector}
    \item{\code{Mitoxantrone_1810}}{a numeric vector}
    \item{\code{Dactinomycin_1811}}{a numeric vector}
    \item{\code{Fludarabine_1813}}{a numeric vector}
    \item{\code{Nelarabine_1814}}{a numeric vector}
    \item{\code{Fulvestrant_1816}}{a numeric vector}
    \item{\code{Vincristine_1818}}{a numeric vector}
    \item{\code{Docetaxel_1819}}{a numeric vector}
    \item{\samp{Podophyllotoxin bromide_1825}}{a numeric vector}
    \item{\code{Dihydrorotenone_1827}}{a numeric vector}
    \item{\code{Gallibiscoquinazole_1830}}{a numeric vector}
    \item{\code{Elephantin_1835}}{a numeric vector}
    \item{\code{Sinularin_1838}}{a numeric vector}
    \item{\code{Sabutoclax_1849}}{a numeric vector}
    \item{\code{LY2109761_1852}}{a numeric vector}
    \item{\samp{OF-1_1853}}{a numeric vector}
    \item{\samp{MN-64_1854}}{a numeric vector}
    \item{\samp{KRAS (G12C) Inhibitor-12_1855}}{a numeric vector}
    \item{\samp{MG-132_1862}}{a numeric vector}
    \item{\samp{BDP-00009066_1866}}{a numeric vector}
    \item{\code{Buparlisib_1873}}{a numeric vector}
    \item{\code{Ulixertinib_1908}}{a numeric vector}
    \item{\code{Venetoclax_1909}}{a numeric vector}
    \item{\code{ABT737_1910}}{a numeric vector}
    \item{\code{Dactinomycin_1911}}{a numeric vector}
    \item{\code{Afuresertib_1912}}{a numeric vector}
    \item{\samp{AGI-5198_1913}}{a numeric vector}
    \item{\code{AZD3759_1915}}{a numeric vector}
    \item{\code{AZD5363_1916}}{a numeric vector}
    \item{\code{AZD6738_1917}}{a numeric vector}
    \item{\code{AZD8186_1918}}{a numeric vector}
    \item{\code{Osimertinib_1919}}{a numeric vector}
    \item{\code{Cediranib_1922}}{a numeric vector}
    \item{\code{Ipatasertib_1924}}{a numeric vector}
    \item{\code{GDC0810_1925}}{a numeric vector}
    \item{\samp{GNE-317_1926}}{a numeric vector}
    \item{\code{GSK2578215A_1927}}{a numeric vector}
    \item{\samp{I-BRD9_1928}}{a numeric vector}
    \item{\samp{Telomerase Inhibitor IX_1930}}{a numeric vector}
    \item{\samp{MIRA-1_1931}}{a numeric vector}
    \item{\samp{NVP-ADW742_1932}}{a numeric vector}
    \item{\code{P22077_1933}}{a numeric vector}
    \item{\code{Savolitinib_1936}}{a numeric vector}
    \item{\samp{UMI-77_1939}}{a numeric vector}
    \item{\code{WIKI4_1940}}{a numeric vector}
    \item{\samp{Sepantronium bromide_1941}}{a numeric vector}
    \item{\code{MIM1_1996}}{a numeric vector}
    \item{\samp{WEHI-539_1997}}{a numeric vector}
    \item{\samp{BPD-00008900_1998}}{a numeric vector}
    \item{\code{Foretinib_2040}}{a numeric vector}
    \item{\samp{BIBR-1532_2043}}{a numeric vector}
    \item{\code{Pyridostatin_2044}}{a numeric vector}
    \item{\samp{AMG-319_2045}}{a numeric vector}
    \item{\samp{MK-8776_2046}}{a numeric vector}
    \item{\code{Ulixertinib_2047}}{a numeric vector}
    \item{\code{Vinorelbine_2048}}{a numeric vector}
    \item{\samp{VX-11e_2096}}{a numeric vector}
    \item{\code{Uprosertib_2106}}{a numeric vector}
    \item{\code{LJI308_2107}}{a numeric vector}
    \item{\code{AZ6102_2109}}{a numeric vector}
    \item{\code{GSK591_2110}}{a numeric vector}
    \item{\code{VE821_2111}}{a numeric vector}
    \item{\code{AZD6482_2169}}{a numeric vector}
    \item{\code{AT13148_2170}}{a numeric vector}
    \item{\samp{BMS-754807_2171}}{a numeric vector}
    \item{\code{JQ1_2172}}{a numeric vector}
  }
}
\references{
Maeser, D., Gruener, R. F., & Huang, R. S. (2021). oncoPredict: an R package for predicting in vivo or cancer patient drug response and biomarkers from cell line screening data. Briefings in bioinformatics, 22(6), bbab260. https://doi.org/10.1093/bib/bbab260
}
\examples{
data(GCAS_drug)
str(GCAS_drug)
}
>>>>>>> 6f4eea48388fdad44a0cc783c245e2ab56dbb4df
