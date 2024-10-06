<<<<<<< HEAD
\name{subtype}
\alias{subtype}
\docType{data}
\title{
Subtypes of al cancer types
}
\description{
Subtype list of al cancer types
}
\usage{data("subtype")}
\format{
  The format is:
List of 21
 $ Bladder cancer      :List of 1
  ..$ BLCA: chr "BLCA"
 $ Breast cancer       :List of 1
  ..$ BRCA:List of 3
  .. ..$ DCIS: chr "DCIS"
  .. ..$ IDC : chr "IDC"
  .. ..$ TNBC: chr "TNBC"
 $ Cervical cancer     :List of 2
  ..$ Premalignant:List of 1
  .. ..$ HSIL: chr "HSIL"
  ..$ CEC         :List of 2
  .. ..$ CIN           : chr "CIN"
  .. ..$ CEC_metastatic: chr "CEC_metastatic"
 $ Colorectal cancer   :List of 2
  ..$ Premalignant:List of 2
  .. ..$ Adenoma: chr "Adenoma"
  .. ..$ polyp  : chr "polyp"
  ..$ CRC         :List of 3
  .. ..$ REC      : chr "REC"
  .. ..$ CRC_liver: chr "CRC_liver"
  .. ..$ CRC_lung : chr "CRC_lung"
 $ Endometrial cancer  :List of 1
  ..$ EMC:List of 2
  .. ..$ EEMC: chr "EEMC"
  .. ..$ SEMC: chr "SEMC"
 $ Esophageal cancer   :List of 1
  ..$ ESC: chr "ESC"
 $ Gastric cancer      :List of 1
  ..$ GAC: chr "GAC"
 $ Glioma              :List of 1
  ..$ Glioma:List of 13
  .. ..$ AA   : chr "AA"
  .. ..$ AO   : chr "AO"
  .. ..$ AOA  : chr "AOA"
  .. ..$ AST  : chr "AST"
  .. ..$ AT/RT: chr "AT/RT"
  .. ..$ EP   : chr "EP"
  .. ..$ MB   : chr "MB"
  .. ..$ OA   : chr "OA"
  .. ..$ OD   : chr "OD"
  .. ..$ PA   : chr "PA"
  .. ..$ PNET : chr "PNET"
  .. ..$ GBM  : chr "GBM"
  .. ..$ rGBM : chr "rGBM"
 $ Head and neck cancer:List of 1
  ..$ HNC:List of 3
  .. ..$ OSCC : chr "OSCC"
  .. ..$ TSCC : chr "TSCC"
  .. ..$ OTSCC: chr "OTSCC"
 $ Kidney cancer       :List of 1
  ..$ RCC: chr "RCC"
 $ Liver cancer        :List of 2
  ..$ Cirrhosis: chr "Cirrhosis"
  ..$ HCC      :List of 2
  .. ..$ CCA         : chr "CCA"
  .. ..$ CirrhosisHCC: chr "CirrhosisHCC"
 $ Lung cancer         :List of 1
  ..$ LC:List of 3
  .. ..$ NSCLC:List of 4
  .. .. ..$ LCC : chr "LCC"
  .. .. ..$ LUAD: chr "LUAD"
  .. .. ..$ LUSC: chr "LUSC"
  .. .. ..$ BAS : chr "BAS"
  .. ..$ CARCI: chr "CARCI"
  .. ..$ LCNE : chr "LCNE"
 $ Skin cancer         :List of 1
  ..$ SC:List of 4
  .. ..$ BCSC         : chr "BCSC"
  .. ..$ MM           : chr "MM"
  .. ..$ MM_metastatic: chr "MM_metastatic"
  .. ..$ SCSC         : chr "SCSC"
 $ Ovarian cancer      :List of 1
  ..$ OC:List of 4
  .. ..$ OCC: chr "OCC"
  .. ..$ MOC: chr "MOC"
  .. ..$ OEC: chr "OEC"
  .. ..$ SOC: chr "SOC"
 $ Pancreatic cancer   :List of 1
  ..$ PANC:List of 2
  .. ..$ PANC_metastaticet: chr "PANC_metastaticet"
  .. ..$ PDAC             : chr "PDAC"
 $ Prostate cancer     :List of 1
  ..$ PCa:List of 2
  .. ..$ PC_metastaticet: chr "PC_metastaticet"
  .. ..$ CRPC           : chr "CRPC"
 $ Thyroid cancer      :List of 1
  ..$ TC:List of 7
  .. ..$ ATC         : chr "ATC"
  .. ..$ FTC         : chr "FTC"
  .. ..$ MTC         : chr "MTC"
  .. ..$ OTC         : chr "OTC"
  .. ..$ PTC         : chr "PTC"
  .. ..$ TC_invasive : chr "TC_invasive"
  .. ..$ PTC_invasive: chr "PTC_invasive"
 $ Osteosarcoma        :List of 1
  ..$ OS: chr "OS"
 $ Myeloma             :List of 2
  ..$ Premalignant:List of 2
  .. ..$ MGUS: chr "MGUS"
  .. ..$ SMM : chr "SMM"
  ..$ MM          :List of 2
  .. ..$ MM         : chr "MM"
  .. ..$ Relapsed MM: chr "Relapsed MM"
 $ Leukemia            :List of 1
  ..$ Leukemia:List of 8
  .. ..$ CLL: chr "CLL"
  .. ..$ PCL: chr "PCL"
  .. ..$ AML: chr "AML"
  .. ..$ ALL:List of 3
  .. .. ..$ ALL  : chr "ALL"
  .. .. ..$ B-ALL: chr "B-ALL"
  .. .. ..$ T-ALL: chr "T-ALL"
  .. ..$ CML: chr "CML"
  .. ..$ MDS: chr "MDS"
  .. ..$ ALT: chr "ALT"
  .. ..$ HCL: chr "HCL"
 $ Lymphoma            :List of 1
  ..$ Lymphoma:List of 15
  .. ..$ WM-BL: chr "WM-BL"
  .. ..$ WM-PC: chr "WM-PC"
  .. ..$ DLBCL: chr "DLBCL"
  .. ..$ FL   : chr "FL"
  .. ..$ BL   : chr "BL"
  .. ..$ cHL  : chr "cHL"
  .. ..$ NLPHL: chr "NLPHL"
  .. ..$ TCRBL: chr "TCRBL"
  .. ..$ CNSI : chr "CNSI"
  .. ..$ MCL  : chr "MCL"
  .. ..$ MZL  :List of 2
  .. .. ..$ MZL : chr "MZL"
  .. .. ..$ NMZL: chr "NMZL"
  .. ..$ BLL  : chr "BLL"
  .. ..$ PEL  : chr "PEL"
  .. ..$ SLVL : chr "SLVL"
  .. ..$ PTCL : chr "PTCL"
}
\examples{
data(subtype)
str(subtype)
}
=======
\name{subtype}
\alias{subtype}
\docType{data}
\title{
Subtypes of al cancer types
}
\description{
Subtype list of al cancer types
}
\usage{data("subtype")}
\format{
  The format is:
List of 21
 $ Bladder cancer      :List of 1
  ..$ BLCA: chr "BLCA"
 $ Breast cancer       :List of 1
  ..$ BRCA:List of 3
  .. ..$ DCIS: chr "DCIS"
  .. ..$ IDC : chr "IDC"
  .. ..$ TNBC: chr "TNBC"
 $ Cervical cancer     :List of 2
  ..$ Premalignant:List of 1
  .. ..$ HSIL: chr "HSIL"
  ..$ CEC         :List of 2
  .. ..$ CIN           : chr "CIN"
  .. ..$ CEC_metastatic: chr "CEC_metastatic"
 $ Colorectal cancer   :List of 2
  ..$ Premalignant:List of 2
  .. ..$ Adenoma: chr "Adenoma"
  .. ..$ polyp  : chr "polyp"
  ..$ CRC         :List of 3
  .. ..$ REC      : chr "REC"
  .. ..$ CRC_liver: chr "CRC_liver"
  .. ..$ CRC_lung : chr "CRC_lung"
 $ Endometrial cancer  :List of 1
  ..$ EMC:List of 2
  .. ..$ EEMC: chr "EEMC"
  .. ..$ SEMC: chr "SEMC"
 $ Esophageal cancer   :List of 1
  ..$ ESC: chr "ESC"
 $ Gastric cancer      :List of 1
  ..$ GAC: chr "GAC"
 $ Glioma              :List of 1
  ..$ Glioma:List of 13
  .. ..$ AA   : chr "AA"
  .. ..$ AO   : chr "AO"
  .. ..$ AOA  : chr "AOA"
  .. ..$ AST  : chr "AST"
  .. ..$ AT/RT: chr "AT/RT"
  .. ..$ EP   : chr "EP"
  .. ..$ MB   : chr "MB"
  .. ..$ OA   : chr "OA"
  .. ..$ OD   : chr "OD"
  .. ..$ PA   : chr "PA"
  .. ..$ PNET : chr "PNET"
  .. ..$ GBM  : chr "GBM"
  .. ..$ rGBM : chr "rGBM"
 $ Head and neck cancer:List of 1
  ..$ HNC:List of 3
  .. ..$ OSCC : chr "OSCC"
  .. ..$ TSCC : chr "TSCC"
  .. ..$ OTSCC: chr "OTSCC"
 $ Kidney cancer       :List of 1
  ..$ RCC: chr "RCC"
 $ Liver cancer        :List of 2
  ..$ Cirrhosis: chr "Cirrhosis"
  ..$ HCC      :List of 2
  .. ..$ CCA         : chr "CCA"
  .. ..$ CirrhosisHCC: chr "CirrhosisHCC"
 $ Lung cancer         :List of 1
  ..$ LC:List of 3
  .. ..$ NSCLC:List of 4
  .. .. ..$ LCC : chr "LCC"
  .. .. ..$ LUAD: chr "LUAD"
  .. .. ..$ LUSC: chr "LUSC"
  .. .. ..$ BAS : chr "BAS"
  .. ..$ CARCI: chr "CARCI"
  .. ..$ LCNE : chr "LCNE"
 $ Skin cancer         :List of 1
  ..$ SC:List of 4
  .. ..$ BCSC         : chr "BCSC"
  .. ..$ MM           : chr "MM"
  .. ..$ MM_metastatic: chr "MM_metastatic"
  .. ..$ SCSC         : chr "SCSC"
 $ Ovarian cancer      :List of 1
  ..$ OC:List of 4
  .. ..$ OCC: chr "OCC"
  .. ..$ MOC: chr "MOC"
  .. ..$ OEC: chr "OEC"
  .. ..$ SOC: chr "SOC"
 $ Pancreatic cancer   :List of 1
  ..$ PANC:List of 2
  .. ..$ PANC_metastaticet: chr "PANC_metastaticet"
  .. ..$ PDAC             : chr "PDAC"
 $ Prostate cancer     :List of 1
  ..$ PCa:List of 2
  .. ..$ PC_metastaticet: chr "PC_metastaticet"
  .. ..$ CRPC           : chr "CRPC"
 $ Thyroid cancer      :List of 1
  ..$ TC:List of 7
  .. ..$ ATC         : chr "ATC"
  .. ..$ FTC         : chr "FTC"
  .. ..$ MTC         : chr "MTC"
  .. ..$ OTC         : chr "OTC"
  .. ..$ PTC         : chr "PTC"
  .. ..$ TC_invasive : chr "TC_invasive"
  .. ..$ PTC_invasive: chr "PTC_invasive"
 $ Osteosarcoma        :List of 1
  ..$ OS: chr "OS"
 $ Myeloma             :List of 2
  ..$ Premalignant:List of 2
  .. ..$ MGUS: chr "MGUS"
  .. ..$ SMM : chr "SMM"
  ..$ MM          :List of 2
  .. ..$ MM         : chr "MM"
  .. ..$ Relapsed MM: chr "Relapsed MM"
 $ Leukemia            :List of 1
  ..$ Leukemia:List of 8
  .. ..$ CLL: chr "CLL"
  .. ..$ PCL: chr "PCL"
  .. ..$ AML: chr "AML"
  .. ..$ ALL:List of 3
  .. .. ..$ ALL  : chr "ALL"
  .. .. ..$ B-ALL: chr "B-ALL"
  .. .. ..$ T-ALL: chr "T-ALL"
  .. ..$ CML: chr "CML"
  .. ..$ MDS: chr "MDS"
  .. ..$ ALT: chr "ALT"
  .. ..$ HCL: chr "HCL"
 $ Lymphoma            :List of 1
  ..$ Lymphoma:List of 15
  .. ..$ WM-BL: chr "WM-BL"
  .. ..$ WM-PC: chr "WM-PC"
  .. ..$ DLBCL: chr "DLBCL"
  .. ..$ FL   : chr "FL"
  .. ..$ BL   : chr "BL"
  .. ..$ cHL  : chr "cHL"
  .. ..$ NLPHL: chr "NLPHL"
  .. ..$ TCRBL: chr "TCRBL"
  .. ..$ CNSI : chr "CNSI"
  .. ..$ MCL  : chr "MCL"
  .. ..$ MZL  :List of 2
  .. .. ..$ MZL : chr "MZL"
  .. .. ..$ NMZL: chr "NMZL"
  .. ..$ BLL  : chr "BLL"
  .. ..$ PEL  : chr "PEL"
  .. ..$ SLVL : chr "SLVL"
  .. ..$ PTCL : chr "PTCL"
}
\examples{
data(subtype)
str(subtype)
}
>>>>>>> 6f4eea48388fdad44a0cc783c245e2ab56dbb4df
