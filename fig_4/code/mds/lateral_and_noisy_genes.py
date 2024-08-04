FIFTY_STRONGEST_PLATELET_GENES = [
    # compiled from pooling the alleged platelets, taking highest 50 genes (50 chose arbitrarily. lowest log_norm_expression (HMGB1) was -8.5)
    'PPBP', 'TMSB4X', 'B2M', 'PF4',
    'FTH1', # https://www.proteinatlas.org/ENSG00000167996-FTH1: "Ferritin heavy chain 1" and "Stores iron in a soluble, non-toxic, readily available form". (https://pubmed.ncbi.nlm.nih.gov/33844865/ seems relevant if curious, but not important for us)
    'RGS18', 'ACTB', 'H3F3A', 'GNG11', 'TUBB1', 'RAP1B', 'SH3BGRL3',
    'TLN1', 
    'CAVIN2',
    'OAZ1', 
    'FTL', # https://www.proteinatlas.org/ENSG00000087086-FTL: "Ferritin light chain" and "Stores iron in a soluble, non-toxic, readily available form".
    'ITM2B', 
    'PTMA', 
    'MYL6', 
    'HIST1H2AC',
    'MYL12A', 
    'RGS10',
    'PRKAR2B', 
    'LIMS1', 
    'CLU', 'NRGN',
    'TSC22D1', 
    'HLA-B',
    'SAT1', 
    'KIF2A', # https://www.proteinatlas.org/ENSG00000068796-KIF2A: "Kinesin family member 2A" and "Plus end-directed microtubule-dependent motor required for normal brain development. May regulate microtubule dynamics during axonal growth. Required for normal progression through mitosis. Required for normal congress of chromosomes at the metaphase plate. Required for normal spindle dynamics during mitosis. Promotes spindle turnover. Implicated in formation of bipolar mitotic spindles. Has microtubule depolymerization activity"
    'OST4', 
    'CCL5', 
    'PLEK', # https://www.proteinatlas.org/ENSG00000115956-PLEK: "Pleckstrin" and "Major protein kinase C substrate of platelets".
    'TAGLN2', 
    'GP9',
    'GSTO1', 
    'RUFY1', 
    'TUBA4A', 'ACRBP',
    'F13A1', # https://www.proteinatlas.org/ENSG00000124491-F13A1: "Coagulation factor XIII A chain".
    'YWHAZ', 
    'TPM3', 
    'CLEC1B',
    'SERF2', 
    'TPM4', 
    'YWHAH', 
    'NCOA4', 
    'SRGN', 
    'HMGB1', 
]
STRONGEST_PLATELET_GENES = [
    *FIFTY_STRONGEST_PLATELET_GENES,

    'TMEM40',
    'HLA-E',
    'PTCRA',
    'VCL', # https://www.proteinatlas.org/ENSG00000035403-VCL: "Vinculin" and "Actin filament (F-actin)-binding protein involved in cell-matrix adhesion and cell-cell adhesion. Regulates cell-surface E-cadherin expression and potentiates mechanosensing by the E-cadherin complex. May also play important roles in cell morphology and locomotion".
    'NT5C3A', # https://www.proteinatlas.org/ENSG00000122643-NT5C3A: "5'-nucleotidase, cytosolic IIIA" and "Nucleotidase which shows specific activity towards cytidine monophosphate (CMP) and 7-methylguanosine monophosphate (m(7)GMP) 1. CMP seems to be the preferred substrate".
    'CALM1', # https://www.proteinatlas.org/ENSG00000198668-CALM1: "Calmodulin 1" and "Calmodulin mediates the control of a large number of enzymes, ion channels, aquaporins and other proteins through calcium-binding. Among the enzymes to be stimulated by the calmodulin-calcium complex are a number of protein kinases and phosphatases. Together with CCP110 and centrin, is involved in a genetic pathway that regulates the centrosome cycle and progression through cytokinesis"
    'H3F3B',
    'NEXN', # https://www.proteinatlas.org/ENSG00000162614-NEXN: "Nexilin F-actin binding protein" and "Involved in regulating cell migration through association with the actin cytoskeleton. Has an essential role in the maintenance of Z line and sarcomere integrity".
    'EIF1', # https://www.proteinatlas.org/ENSG00000173812-EIF1: "Eukaryotic translation initiation factor 1" and "Necessary for scanning and involved in initiation site selection. Promotes the assembly of 48S ribosomal complexes at the authentic initiation codon of a conventional capped mRNA".
    'C2orf88', # https://www.proteinatlas.org/ENSG00000187699-C2orf88: "Chromosome 2 open reading frame 88" and "Binds to type I regulatory subunits of protein kinase A (PKA-RI) and may anchor/target them to the plasma membrane"
    'HLA-A',
    'MTPN', # https://www.proteinatlas.org/ENSG00000105887-MTPN: "Myotrophin" and "Promotes dimerization of NF-kappa-B subunits and regulates NF-kappa-B transcription factor activity (By similarity). Plays a role in the regulation of the growth of actin filaments. Inhibits the activity of the F-actin-capping protein complex formed by the CAPZA1 and CAPZB heterodimer".
    'MPP1', # https://www.proteinatlas.org/ENSG00000130830-MPP1: "Membrane palmitoylated protein 1" and "Essential regulator of neutrophil polarity. Regulates neutrophil polarization by regulating AKT1 phosphorylation through a mechanism that is independent of PIK3CG activity" and "This gene encodes the prototype of the membrane-associated guanylate kinase (MAGUK) family proteins. MAGUKs interact with the cytoskeleton and regulate cell proliferation, signaling pathways, and intercellular junctions. The encoded protein is an extensively palmitoylated membrane phosphoprotein".
    'ARPC1B', # https://www.proteinatlas.org/ENSG00000130429-ARPC1B: "Actin related protein 2/3 complex subunit 1B" and "Component of the Arp2/3 complex, a multiprotein complex that mediates actin polymerization upon stimulation by nucleation-promoting factor (NPF) 1, 2. The Arp2/3 complex mediates the formation of branched actin networks in the cytoplasm, providing the force for cell motility 3, 4. In addition to its role in the cytoplasmic cytoskeleton, the Arp2/3 complex also promotes actin polymerization in the nucleus, thereby regulating gene transcription and repair of damaged DNA 5. The Arp2/3 complex promotes homologous recombination (HR) repair in response to DNA damage by promoting nuclear actin polymerization, leading to drive motility of double-strand breaks (DSBs)" and "This gene encodes one of seven subunits of the human Arp2/3 protein complex".
    'ACTG1',
    'ARPC5', # https://www.proteinatlas.org/ENSG00000162704-ARPC5: "Actin related protein 2/3 complex subunit 5" and "This gene encodes one of seven subunits of the human Arp2/3 protein complex"
    'MAX', # https://www.proteinatlas.org/ENSG00000125952-MAX: "MYC associated factor X" and "Transcription regulator. Forms a sequence-specific DNA-binding protein complex with MYC or MAD which recognizes the core sequence 5'-CAC[GA]TG-3'. The MYC:MAX complex is a transcriptional activator, whereas the MAD:MAX complex is a repressor. May repress transcription via the recruitment of a chromatin remodeling complex containing H3 'Lys-9' histone methyltransferase activity. Represses MYC transcriptional activity from E-box elements".
    'DAB2', # https://www.proteinatlas.org/ENSG00000153071-DAB2: "DAB adaptor protein 2" and "Adapter protein that functions as clathrin-associated sorting protein (CLASP) required for clathrin-mediated endocytosis of selected cargo proteins. Can bind and assemble clathrin, and binds simultaneously to phosphatidylinositol 4,5-bisphosphate (PtdIns(4,5)P2) and cargos containing non-phosphorylated NPXY internalization motifs, such as the LDL receptor, to recruit them to clathrin-coated pits."
    'TIMP1', # https://www.proteinatlas.org/ENSG00000102265-TIMP1: "TIMP metallopeptidase inhibitor 1" and "Metalloproteinase inhibitor that functions by forming one to one complexes with target metalloproteinases, such as collagenases, and irreversibly inactivates them by binding to their catalytic zinc cofactor".
    'CD9', # https://www.proteinatlas.org/ENSG00000010278-CD9: "Integral membrane protein associated with integrins, which regulates different processes, such as sperm-egg fusion, platelet activation and aggregation, and cell adhesion".
    'UQCRH', # https://www.proteinatlas.org/ENSG00000173660-UQCRH: "Ubiquinol-cytochrome c reductase hinge protein" and "Component of the ubiquinol-cytochrome c oxidoreductase, a multisubunit transmembrane complex that is part of the mitochondrial electron transport chain which drives oxidative phosphorylation"
    'RIOK3', # https://www.proteinatlas.org/ENSG00000101782-RIOK3: "RIO kinase 3" and "Involved in regulation of type I interferon (IFN)-dependent immune response which plays a critical role in the innate immune response against DNA and RNA viruses".
    'MMD',
    'ODC1', # https://www.proteinatlas.org/ENSG00000115758-ODC1: "Ornithine decarboxylase 1" and "Catalyzes the first and rate-limiting step of polyamine biosynthesis that converts ornithine into putrescine, which is the precursor for the polyamines, spermidine and spermine. Polyamines are essential for cell proliferation and are implicated in cellular processes, ranging from DNA replication to apoptosis".
    'TPM1', # https://www.proteinatlas.org/ENSG00000140416-TPM1: "Tropomyosin 1" and "Binds to actin filaments in muscle and non-muscle cells 1. Plays a central role, in association with the troponin complex, in the calcium dependent regulation of vertebrate striated muscle contraction 2. Smooth muscle contraction is regulated by interaction with caldesmon. In non-muscle cells is implicated in stabilizing cytoskeleton actin filaments".
    'GAPDH', # https://www.proteinatlas.org/ENSG00000111640-GAPDH: "Glyceraldehyde-3-phosphate dehydrogenase" and "Has both glyceraldehyde-3-phosphate dehydrogenase and nitrosylase activities, thereby playing a role in glycolysis and nuclear functions, respectively".
    'SH3BP5', # https://www.proteinatlas.org/ENSG00000131370-SH3BP5: "SH3 domain binding protein 5" and "Functions as guanine nucleotide exchange factor (GEF) with specificity for RAB11A and RAB25 1, 2. Inhibits the auto- and transphosphorylation activity of BTK. Plays a negative regulatory role in BTK-related cytoplasmic signaling in B-cells. May be involved in BCR-induced apoptotic cell death"
    'CDKN2D', # https://www.proteinatlas.org/ENSG00000129355-CDKN2D: "Cyclin dependent kinase inhibitor 2D" and "Interacts strongly with CDK4 and CDK6 and inhibits them".
    'SEPT7', # https://www.proteinatlas.org/ENSG00000122545-SEPTIN7: "Septin 7" and "Filament-forming cytoskeletal GTPase. Required for normal organization of the actin cytoskeleton. Required for normal progress through mitosis. Involved in cytokinesis. Required for normal association of CENPE with the kinetochore. Plays a role in ciliogenesis and collective cell movements".
    'TREML1',
    'LGALSL',
    'ARHGAP18', # https://www.proteinatlas.org/ENSG00000146376-ARHGAP18: "Rho GTPase activating protein 18" and "Rho GTPase activating protein that suppresses F-actin polymerization by inhibiting Rho. Rho GTPase activating proteins act by converting Rho-type GTPases to an inactive GDP-bound state 1. Plays a key role in tissue tension and 3D tissue shape by regulating cortical actomyosin network formation"
    'PRDX6', # https://www.proteinatlas.org/ENSG00000117592-PRDX6: "Peroxiredoxin 6" and "Thiol-specific peroxidase that catalyzes the reduction of hydrogen peroxide and organic hydroperoxides to water and alcohols, respectively"
    'CALM3', # https://www.proteinatlas.org/ENSG00000160014-CALM3: "Calmodulin 3" and "Calmodulin mediates the control of a large number of enzymes, ion channels, aquaporins and other proteins through calcium-binding. Is a regulator of voltage-dependent L-type calcium channels 1. Among the enzymes to be stimulated by the calmodulin-calcium complex are a number of protein kinases and phosphatases. Together with CCP110 and centrin, is involved in a genetic pathway that regulates the centrosome cycle and progression through cytokinesis".
    'RAB32', # https://www.proteinatlas.org/ENSG00000118508-RAB32: "RAB32, member RAS oncogene family" and "Acts as an A-kinase anchoring protein by binding to the type II regulatory subunit of protein kinase A and anchoring it to the mitochondrion. Also involved in synchronization of mitochondrial fission 1. Plays a role in the maturation of phagosomes that engulf pathogens, such as S.aureus and M.tuberculosis".
    'ARL6IP5', # https://www.proteinatlas.org/ENSG00000144746-ARL6IP5: "ADP ribosylation factor like GTPase 6 interacting protein 5" and "Regulates intracellular concentrations of taurine and glutamate. Negatively modulates SLC1A1/EAAC1 glutamate transport activity by decreasing its affinity for glutamate in a PKC activity-dependent manner. Plays a role in the retention of SLC1A1/EAAC1 in the endoplasmic reticulum".
    'DYNLL1', # https://www.proteinatlas.org/ENSG00000088986-DYNLL1: "Dynein light chain LC8-type 1" and "Acts as one of several non-catalytic accessory components of the cytoplasmic dynein 1 complex that are thought to be involved in linking dynein to cargos and to adapter proteins that regulate dynein function. Cytoplasmic dynein 1 acts as a motor for the intracellular retrograde motility of vesicles and organelles along microtubules. May play a role in changing or maintaining the spatial distribution of cytoskeletal structures"
    'ATP5F1E', 
    'MTURN', # https://www.proteinatlas.org/ENSG00000180354-MTURN: "Maturin, neural progenitor differentiation regulator homolog" and "Promotes megakaryocyte differentiation by enhancing ERK and JNK signaling as well as up-regulating RUNX1 and FLI1 expression".
    'RBX1', # https://www.proteinatlas.org/ENSG00000100387-RBX1: "Ring-box 1" and "E3 ubiquitin ligase component of multiple cullin-RING-based E3 ubiquitin-protein ligase (CRLs) complexes which mediate the ubiquitination and subsequent proteasomal degradation of target proteins, including proteins involved in cell cycle progression, signal transduction, transcription and transcription-coupled nucleotide excision repair".
    'GNAS', # https://www.proteinatlas.org/ENSG00000087460-GNAS: "GNAS complex locus" and "Guanine nucleotide-binding proteins (G proteins) function as transducers in numerous signaling pathways controlled by G protein-coupled receptors (GPCRs) 1. Signaling involves the activation of adenylyl cyclases, resulting in increased levels of the signaling molecule cAMP 2, 3. GNAS functions downstream of several GPCRs, including beta-adrenergic receptors 4. Stimulates the Ras signaling pathway via RAPGEF2" and "This locus has a highly complex imprinted expression pattern. It gives rise to maternally, paternally, and biallelically expressed transcripts that are derived from four alternative promoters and 5' exons".
    'CFL1', # https://www.proteinatlas.org/ENSG00000172757-CFL1: "Cofilin 1" and "Binds to F-actin and exhibits pH-sensitive F-actin depolymerizing activity".
    'RAB11A', # https://www.proteinatlas.org/ENSG00000103769-RAB11A: "RAB11A, member RAS oncogene family" and "The small GTPases Rab are key regulators of intracellular membrane trafficking, from the formation of transport vesicles to their fusion with membranes. Rabs cycle between an inactive GDP-bound form and an active GTP-bound form that is able to recruit to membranes different set of downstream effectors directly responsible for vesicle formation, movement, tethering and fusion. The small Rab GTPase RAB11A regulates endocytic recycling. Acts as a major regulator of membrane delivery during cytokinesis."
    'PPDPF', # https://www.proteinatlas.org/ENSG00000125534-PPDPF: "Pancreatic progenitor cell differentiation and proliferation factor" and "Probable regulator of exocrine pancreas development".
    'RAB27B', # https://www.proteinatlas.org/ENSG00000041353-RAB27B: "RAB27B, member RAS oncogene family" and "Small GTPase which cycles between active GTP-bound and inactive GDP-bound states. In its active state, binds to a variety of effector proteins to regulate homeostasis of late endocytic pathway, including endosomal positioning, maturation and secretion"
    'ARHGDIB', # https://www.proteinatlas.org/ENSG00000111348-ARHGDIB: "Rho GDP dissociation inhibitor beta" and "Regulates the GDP/GTP exchange reaction of the Rho proteins by inhibiting the dissociation of GDP from them, and the subsequent binding of GTP to them 1, 2. Regulates reorganization of the actin cytoskeleton mediated by Rho family members".
    'CTSA', # https://www.proteinatlas.org/ENSG00000064601-CTSA: "Cathepsin A" and "Protective protein appears to be essential for both the activity of beta-galactosidase and neuraminidase, it associates with these enzymes and exerts a protective function necessary for their stability and activity. This protein is also a carboxypeptidase and can deamidate tachykinins".
    'EIF2AK1',
    'PLEKHO1', # https://www.proteinatlas.org/ENSG00000023902-PLEKHO1: "Pleckstrin homology domain containing O1" and "Plays a role in the regulation of the actin cytoskeleton through its interactions with actin capping protein (CP). May function to target CK2 to the plasma membrane thereby serving as an adapter to facilitate the phosphorylation of CP by protein kinase 2 (CK2). Appears to target ATM to the plasma membrane. Appears to also inhibit tumor cell growth by inhibiting AKT-mediated cell-survival. Also implicated in PI3K-regulated muscle differentiation, the regulation of AP-1 activity (plasma membrane bound AP-1 regulator that translocates to the nucleus) and the promotion of apoptosis induced by tumor necrosis factor TNF. When bound to PKB, it inhibits it probably by decreasing PKB level of phosphorylation".
]


PROBABLY_PLATELET_SPECIFIC_GENES = [
    'PF4', 'PPBP', # these are assumptions.

    # didn't test this quantitatively, but comparison of my alleged platelet transcriptional profile to https://ashpublications.org/blood/article/118/14/e101/28765/Genome-wide-RNA-seq-analysis-of-human-and-mouse passes my sanity check. specifically, it also seems there that the amount of ribosomal protein mRNAs is very low.

    # this list was compiled by pooling downsampled UMIs of ~900 alleged platelets (based on downsampled PF4+PPBP)
    # (platelet_mask = (
    #     (total_platelet_marker_gene_downsampled_umis_of_cells >= 15)
    #     & (clean_ad.obs['log_norm_num_of_cell_ranger_reported_umis'] <= -1)
    #     & (ribosomal_protein_gene_umi_fraction <= 0.05)
    # ))
    # and comparing to probably not platelets (based on PF4+PPBP)
    # (probably_not_platelet_mask = (
    #     (mc.ut.to_numpy_vector(clean_ad.X[:, platelet_marker_gene_indices].sum(axis=1)) == 0)
    #     & (ribosomal_protein_gene_umi_fraction >= 0.1)
    # ))
    

    # 230304: all of the following genes are at least 4 times higher in platelets compared to MKPs.

    'GP9', 'CLEC1B', 'TMEM40',
    # 'PTCRA', # can be not-low in monocytes
    'TREML1', 'LGALSL',
    'CLDN5', # can be high in endothel, but i guess we don't care about endothel that much
    # 'SMIM5', # can be not-low in monocytes
    'CMTM5',
    # 'AP001189.1', # feels unsafe to use it.
    'ENKUR', # kind of can be not-low in BEMP
    'CTTN', # can be high in endothel, but i guess we don't care about endothel that much
    'PDE5A', 'LCN2',
    # 'AC147651.1', # feels unsafe to use it.
    'PDGFA', # can be not-low in endothel, but i guess we don't care about endothel that much
    'TNNC2',
    'C19orf33', # can be not-low in endothel, but i guess we don't care about endothel that much
    # 'AC090409.1', # feels unsafe to use it.
    'ITGB3', # can be not-low in endothel, but i guess we don't care about endothel that much
    'GP1BA', 'HGD', 'TUBA8',
    # 'CDKN1A', # can be not-low in monocytes
    'LGALS12', 'MFAP3L', 'DENND2C',
    'FRMD3', # kind of can be not-low in monocytes
    # 'ANKRD9', # can be not-low in EP
    'SLFN14',
    # 'ABCC3', # can be not-low in monocytes
    'AQP10',
    # 'SCN1B', # can be not-low in CLP
    'PROS1', # can be not-low in endothel, but i guess we don't care about endothel that much
    'TRAPPC3L', 'DMTN',
    # 'AC000093.1', # feels unsafe to use it.
    # 'GNAZ', # can be not-low in B
    'CXCL5',
    # 'TSPAN18', # can be not-low in NKT
]

MAYBE_NOISY_6P_HLA_GENE_NAMES = [
    # this is only one cluster of HLA genes. there are more nearby.
    'HLA-DRA', 'HLA-DRB5', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DQB1-AS1', 'HLA-DQA2', 'HLA-DQB2', 'HLA-DOB',
]
NOISY_6P_HIST_GENE_NAMES = [
    'HIST1H1A', 'HIST1H3A', 'HIST1H4A', 'HIST1H4B',
    'HIST1H3B', 'HIST1H2AB', 'HIST1H2BB', 'HIST1H3C', 'HIST1H1C',
    'HFE', 'HIST1H4C', 'HIST1H1T', 'HIST1H2BC', 'HIST1H2AC',
    'HIST1H1E', 'HIST1H2BD', 'AL353759.1', 'HIST1H2BE', 'HIST1H4D',
    'AL031777.3', 'HIST1H2AD', 'HIST1H2BF', 'HIST1H4E', 'HIST1H2BG',
    'HIST1H2AE', 'HIST1H3E', 'HIST1H1D', 'HIST1H4F', 'HIST1H3F',
    'HIST1H2BH', 'HIST1H3G', 'HIST1H2BI', 'HIST1H4H',
]

IFN_MODULE_IN_HIGH_IFN_METACELLS = [
    # this list was compiled by clustering genes over metacells, but only 'Dissimilar-high-IFN' metacells, and only genes with log_norm_expression_range >= 2.5.
    'IFIT2', 'PARP14', 'APOL6', 'IFI35', 'PLSCR1', 'STAT2',
    'SP110', 'XRN1', 'OAS2', 'MX1', 'SAMD9L', 'MX2',
    'OASL', 'DDX60L', 'SAMD9',
    'HERC5', 'ISG15', 'LY6E',
    'IFI6', 'DDX58', 'IFIT3', 'IFIT1',
    'RSAD2', 'IFIH1', 'OAS3',
    'IFI44L', 'DDX60', 'DTX3L', 'PARP9', 'IFI44', 'XAF1',
]

PRESUMABLY_UNIVERSALLY_HIGH_IN_DISSIMILAR_HIGH_IFN = [
    'IFI6', 'OAS3', 'IFI44', 'STAT1', 'OAS2', 'DDX60', 'PLSCR1', 'SAMD9L', 'XAF1', 'IFIH1', 'EPSTI1', 'PARP14', 'BST2', 'IFITM3', 'ISG15', 'OAS1', 'IFIT3', 'LAP3', 'DDX58', 'MX1', 'IFI44L', 'IFIT1',
]
PRESUMABLY_COMMONLY_HIGH_IN_DISSIMILAR_HIGH_IFN = [
    'HLA-C', 'GBP1', 'PHF11', 'PARP9', 'CTCF', 'XRN1', 'DTX3L', 'LY6E', 'IFIT5', 'DDX60L', 'RNF213', 'HERC5', 'SAMD9', 'IRF7', 'MX2', 'MT2A', 'ADAR', 'EIF2AK2', 'RSAD2', 'IFIT2', 'EIF1AY', 'HIST1H1E', 'EGR1', 'KDM5D', 'CCND2', 'HLA-B', 'SPEN', 'TAP2', 'GPR180', 'IER2',
]
PRESUMABLY_HIGH_IN_DISSIMILAR_HIGH_IFN_AT_LEAST_IN_A_MEBEMP_M_AND_MPP_CLUSTER = [
    'SLFN5', 'FUS', 'PARP12', 'IFI35', 'CHMP5', 'SPATS2L', 'N4BP1', 'TRIM14', 'UBE2D3', 'NMI', 'DRAP1', 'ODF3B', 'SP110', 'GBP2', 'B2M', 'ISG20', 'CMPK2', 'DDX3Y', 'DNAJA1', 'SMCHD1',
]

RIBOSOMAL_PROTEIN_GENE_NAMES = [
    # 'RP[SL][0-9]*$'
    'RPS8',
    'RPS7',
    'RPS9',
    'RPS6',
    'RPS3',
    'RPS2',
    'RPS5',
    'RPL28',
    'RPS13',
    'RPS25',
    'RPS24',
    'RPS26',
    'RPL41',
    'RPL6',
    'RPL21',
    'RPS29',
    'RPL4',
    'RPS17',
    'RPL13',
    'RPL26',
    'RPL23',
    'RPL19',
    'RPL27',
    'RPL38',
    'RPL17',
    'RPS21',
    'RPS15',
    'RPL36',
    'RPS28',
    'RPS16',
    'RPS19',
    'RPL18',
    'RPS11',
    'RPL3',
    'RPS27',
    'RPL22',
    'RPL11',
    'RPL5',
    'RPL31',
    'RPL32',
    'RPL15',
    'RPS12',
    'RPL14',
    'RPL29',
    'RPL24',
    'RPL9',
    'RPL34',
    'RPL37',
    'RPS23',
    'RPS14',
    'RPS18',
    'RPS10',
    'RPL39',
    'RPL10',
    'RPS20',
    'RPL7',
    'RPL30',
    'RPL8',
    'RPL35',
    'RPL12',

    'RPS15A',
    'RPL23A',
    'RPL18A',
    'RPL13A',
    'RPL7A',
    'RPL27A',
    'RPL36A',
    'RPL37A',
    'RPS27A',
    'RPL35A',
    'RPS3A',
    'RPL10A',

    'RPSA', # seems like this is RPS1 (e.g., https://www.nature.com/articles/s41598-019-44013-9: "ribosomal protein S1 (RpsA)")


    # https://www.proteinatlas.org/ENSG00000198034-RPS4X: Ribosomal protein S4 is the only ribosomal protein known to be encoded by more than one gene, namely this gene and ribosomal protein S4, Y-linked (RPS4Y). The 2 isoforms encoded by these genes are not identical, but are functionally equivalent.
    'RPS4X',
    'RPS4Y1',
    'RPS4Y2', # https://www.proteinatlas.org/ENSG00000280969-RPS4Y2: "The protein encoded by this gene is a ribosomal protein that is highly similar to RPS4Y1."

    'RPLP0', # https://www.proteinatlas.org/ENSG00000089157-RPLP0: "Ribosomal protein lateral stalk subunit P0" and "This gene encodes a ribosomal protein that is a component of the 60S subunit."
    'RPLP1', # https://www.proteinatlas.org/ENSG00000137818-RPLP1: "Ribosomal protein lateral stalk subunit P1" and "This gene encodes a ribosomal phosphoprotein that is a component of the 60S subunit."
    'RPLP2', # https://www.proteinatlas.org/ENSG00000177600-RPLP2: "Ribosomal protein lateral stalk subunit P2" and "This gene encodes a ribosomal phosphoprotein that is a component of the 60S subunit."

    'RPL36AL', # https://www.proteinatlas.org/ENSG00000165502-RPL36AL: "Ribosomal protein L36a like" and "This gene encodes a ribosomal protein that is a component of the 60S subunit."


    'RPL22L1', 'RPS27L', 'RPL26L1', 'RPL39L', # decided to not exclude (out of the ribosomal protein gene list) these in the end. see below.
]

REPLICATION_DEPENDENT_HISTONE_GENE_NAMES = [
    # 221210: HIST1H1D looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
    'HIST1H1D', # https://www.proteinatlas.org/ENSG00000124575-H1-3: "H1.3 linker histone, cluster member" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H1 family."

    'HIST2H2AA3', # https://www.proteinatlas.org/ENSG00000203812-H2AC18: "H2A clustered histone 18" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H2A family"

    'HIST1H3G', # https://www.proteinatlas.org/ENSG00000273983-H3C8: "H3 clustered histone 8" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H3 family."

    'HIST1H3H', # https://www.proteinatlas.org/ENSG00000278828-H3C10: "H3 clustered histone 10" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H3 family"

    'HIST1H2BJ', # https://www.proteinatlas.org/ENSG00000124635-H2BC11: "H2B clustered histone 11" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H2B family"

    'HIST2H2BF', # https://www.proteinatlas.org/ENSG00000203814-H2BC18: "H2B clustered histone 18" and "This gene encodes a replication-dependent histone that is a member of the histone H2B family"

    'HIST1H2AG', # https://www.proteinatlas.org/ENSG00000196787-H2AC11: "H2A clustered histone 11" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H2A family"

    'HIST1H2BF', # https://www.proteinatlas.org/ENSG00000277224-H2BC7: "H2B clustered histone 7" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H2B family"

    # 221210: HIST1H4C looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
    'HIST1H4C', # https://www.proteinatlas.org/ENSG00000197061-H4C3: "H4 clustered histone 3" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H4 family"

    'HIST1H4E', # https://www.proteinatlas.org/ENSG00000276966-H4C5: "H4 clustered histone 5" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H4 family"

    # 221210: HIST1H1E looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
    'HIST1H1E', # https://www.proteinatlas.org/ENSG00000168298-H1-4: "H1.4 linker histone, cluster member" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H1 family"

    # 221210: HIST1H2AC looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
    'HIST1H2AC', # https://www.proteinatlas.org/ENSG00000180573-H2AC6: "H2A clustered histone 6" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H2A family". also, in metacells with very high PF4 and/or PPBP, HIST1H2AC is also very high

    # 221210: HIST1H4H looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
    'HIST1H4H', # https://www.proteinatlas.org/ENSG00000158406-H4C8: "H4 clustered histone 8" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H4 family"

    # 221210: HIST1H1B looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
    'HIST1H1B', # https://www.proteinatlas.org/ENSG00000184357-H1-5: "H1.5 linker histone, cluster member" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H1 family."

    'HIST1H2BC', # https://www.proteinatlas.org/ENSG00000180596-H2BC4: "H2B clustered histone 4" and "The main transcript variant of this gene is intronless and encodes a replication-dependent histone that is a member of the histone H2B family"
    'HIST1H1C', # https://www.proteinatlas.org/ENSG00000187837-H1-2: "H1.2 linker histone, cluster member" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H1 family"

    'HIST1H2BD', # https://www.proteinatlas.org/ENSG00000158373-H2BC5: "H2B clustered histone 5" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H2B family".

    'HIST2H2AB', # https://www.proteinatlas.org/ENSG00000184270-H2AC21: "H2A clustered histone 21" and "This gene is intronless and encodes a replication-dependent histone that is a member of the histone H2A family"
]

GENES_EXCLUDED_IN_THE_METACELL_VIGNETTE = [
    # ugh. i don't have good excluses to mark the following as lateral, but i just assume metacells vignette is an expert.
    # the following genes appear in https://github.com/tanaylab/metacells/blob/master/vignettes/Metacells_Vignette.ipynb in the excluded gene list.
    # https://www.proteinatlas.org/ENSG00000132740-IGHMBP2:
    # "5' to 3' helicase that unwinds RNA and DNA duplices in an ATP-dependent reaction. Acts as a transcription regulator." and
    # "This gene encodes a helicase superfamily member that binds a specific DNA sequence from the immunoglobulin mu chain switch region."
    'IGHMBP2',
    # https://www.proteinatlas.org/ENSG00000128322-IGLL1: "The preB cell receptor is found on the surface of proB and preB cells, where
    # it is involved in transduction of signals for cellular proliferation, differentiation from the proB cell to the preB cell stage,
    # allelic exclusion at the Ig heavy chain gene locus, and promotion of Ig light chain gene rearrangements."
    'IGLL1',
    # https://www.proteinatlas.org/ENSG00000254709-IGLL5: "This gene encodes one of the immunoglobulin lambda-like polypeptides.
    # It is located within the immunoglobulin lambda locus but it does not require somatic rearrangement for expression."
    'IGLL5',
    # https://www.genecards.org/cgi-bin/carddisp.pl?gene=IGLON5: "Predicted to be located in extracellular region."
    # https://www.frontiersin.org/articles/10.3389/fimmu.2022.852215/full: "The IgLON family consists of five genes: Lsamp, Ntm, Opcml, Negr1,
    # and Iglon5. It is well known for being involved in the process of neuronal adhesion, neurogenesis and neuroplasticity and is strictly
    # regulated by metalloproteinase activity on the surface of mature cortical neurons"
    'IGLON5',

    # 221210: TMSB4X is in https://github.com/tanaylab/metacells/blob/master/vignettes/Metacells_Vignette.ipynb as one of the well known genes to exclude, and in my data also its levels are very high. also, its presumable function seems relatively lateral. so i guess marking it as lateral makes sense. also, looked batchy.
    'TMSB4X', # https://www.proteinatlas.org/ENSG00000205542-TMSB4X: "Thymosin beta 4 X-linked" and "Plays an important role in the organization of the cytoskeleton 1, 2. Binds to and sequesters actin monomers (G actin) and therefore inhibits actin polymerization" and "This gene encodes an actin sequestering protein which plays a role in regulation of actin polymerization. The protein is also involved in cell proliferation, migration, and differentiation. This gene escapes X inactivation and has a homolog on chromosome Y."
]

GENES_WITH_MORE_THAN_4_ULT_ILL_FOLD = [
    'AC005921.2', 'AC008035.1', 'AC010642.2', 'AC010883.1',
    'AC015849.3', 'AC087190.1', 'AC104532.2', 'ACOT2', 'AL022238.2',
    'AL049840.1', 'AL591767.1', 'APOLD1', 'ARF5', 'ARL5A', 'CALHM6',
    'CRIP1', 'DGCR6L', 'DHRS4L2', 'EEF1G', 'EGLN2', 'EIF3CL', 'F8A1',
    'FBXL6', 'GNG10', 'HIST1H4J', 'HIST2H2AA3', 'HIST2H4B', 'HLA-DRB5',
    'IGLC2', 'IGLC3', 'KCNC1', 'KMT2E-AS1', 'LINC00513', 'LINC02256',
    'MIF', 'MT-ATP8', 'MT-CO2', 'MT-ND4L', 'NDUFA7', 'PAGR1',
    'PDCD4-AS1', 'POLD4', 'RAB24', 'RBAK-RBAKDN', 'RPL17', 'RPL23A',
    'RPL36A', 'RPS17', 'SF3B4', 'TCIRG1', 'TGFB1', 'TMEM191B',
    'TMPRSS9', 'TRIM73', 'TTC4', 'ZP3',
]
GENES_WITH_HIGH_ULT_ILL_LOG_RATIO_STD_AND_AT_LEAST_SOMETIME_ULT_ILL_DIFF_EXP = [
    'ABHD14A-ACY1', 'ABHD16A', 'AC000093.1', 'AC002464.1',
    'AC004540.2', 'AC004556.1', 'AC004771.4', 'AC005670.2',
    'AC005921.2', 'AC006064.4', 'AC007014.1', 'AC007036.1',
    'AC007325.4', 'AC007952.4', 'AC008610.1', 'AC008676.1',
    'AC010761.1', 'AC010880.1', 'AC010883.1', 'AC010931.2',
    'AC010978.1', 'AC012306.2', 'AC012615.2', 'AC012645.2',
    'AC012645.3', 'AC015726.1', 'AC015849.3', 'AC016735.1',
    'AC018695.2', 'AC020656.1', 'AC020910.4', 'AC020911.2',
    'AC023034.1', 'AC023983.1', 'AC025171.3', 'AC026471.3',
    'AC036214.1', 'AC040977.1', 'AC048382.6', 'AC068491.3',
    'AC069277.1', 'AC074183.1', 'AC087190.1', 'AC091057.6',
    'AC092069.1', 'AC092140.1', 'AC092171.3', 'AC092301.1',
    'AC092691.1', 'AC092723.1', 'AC093274.1', 'AC093495.1',
    'AC097532.2', 'AC099811.4', 'AC100793.2', 'AC104041.1',
    'AC104248.1', 'AC104532.2', 'AC107871.1', 'AC109597.2',
    'AC112229.3', 'AC114490.2', 'AC116914.2', 'AC123768.2',
    'AC125603.1', 'AC127521.1', 'AC132192.1', 'AC133919.1',
    'AC136475.5', 'AC140725.1', 'AC146944.4', 'AC211486.2',
    'AC211486.5', 'AC239799.2', 'AC239868.3', 'AC240274.1',
    'AC245595.1', 'AC246785.3', 'AC246793.1', 'ACOT2', 'ACY1',
    'AD000090.1', 'ADSL', 'AF064858.1', 'AF127577.1', 'AGAP9', 'AK1',
    'AKAP2', 'AKT1', 'AL022238.2', 'AL031714.1', 'AL031777.3',
    'AL032821.1', 'AL035078.1', 'AL049840.1', 'AL121944.1',
    'AL133453.1', 'AL133523.1', 'AL135925.1', 'AL136038.5',
    'AL136295.2', 'AL136454.1', 'AL139353.1', 'AL157938.2',
    'AL158154.2', 'AL158212.2', 'AL160191.1', 'AL354733.3',
    'AL355432.1', 'AL355816.2', 'AL357093.1', 'AL359636.1',
    'AL360012.1', 'AL390728.5', 'AL445686.2', 'AL451165.2',
    'AL512770.1', 'AL590764.1', 'AL590867.1', 'AL591767.1',
    'AL603832.1', 'AL662844.4', 'AL731557.1', 'AL929472.2', 'ALDOA',
    'AMIGO3', 'AP000763.3', 'AP001160.1', 'AP001267.2', 'AP002433.1',
    'AP002449.1', 'AP003392.1', 'APOLD1', 'ARAP1', 'ARF4-AS1', 'ARF5',
    'ARPC1B', 'ATN1', 'ATP5F1D', 'ATP5MGL', 'ATP6V1B1-AS1', 'ATXN2L',
    'ATXN8OS', 
    ###'AVP', 
    'BCKDHA', 'BOLA2B', 'BTG2', 'BX284668.5',
    'C11orf98', 'C16orf95', 'C17orf49', 'C18orf32', 'C20orf144',
    'C5orf66', 'CALHM2', 'CALHM6', 'CAMTA1-DT', 'CBWD6', 'CCBE1',
    'CCDC39', 'CCDC71', 'CCZ1B', 'CD24', 'CD70', 'CD8B', 'CDK11A',
    'CDRT4', 'CDV3', 'CEBPA', 'CERKL', 'CERS2', 'CHKB', 'CHMP1B-AS1',
    'CMTM5', 'CNTN2', 'COCH', 'COL6A2', 'COL6A3', 'COMT', 'COX16',
    'CPLX1', 'CPNE7', 'CPT1B', 'CRADD', 'CRIP1', 'CRIP2', 'CRTC2',
    'CRYBB2', 'CSNK2A3', 'CTAG2', 'CXorf40A', 'CXorf40B', 'CYBA',
    'DDIT4-AS1', 'DDR1-DT', 'DDTL', 'DGCR6', 'DHRS4', 'DHRS4L2',
    'DND1', 'DUSP2', 'DYTN', 'EBI3', 'EDDM13', 'EDF1', 'EEF1G',
    'EGLN2', 'EIF2S3B', 'EIF3CL', 'EIF4EBP3', 'EIF5AL1', 'ESCO2',
    'F8A1', 'F8A3', 'FALEC', 'FAM153C', 'FAM157C', 'FAM234A',
    'FAM239A', 'FAM239B', 'FAM72D', 'FBXL6', 'FCGR3B', 'FDX2',
    'FKBP11', 'FOXJ1', 'FP565260.1', 'FRS3', 'FSIP2', 'FXYD1',
    'GABARAP', 'GATD3A', 'GBP7', 'GCC2-AS1', 'GDF10', 'GGT1', 'GHRL',
    'GNG10', 'GNG8', 'GOLGA6L10', 'GOLGA6L9', 'GOLGA8B', 'GOLGA8Q',
    'GOLGA8S', 'GPM6A', 'GPR89A', 'GPX2', 'GRK2', 'GSDMD', 'GTF2IRD1',
    'GTF2IRD2B', 'GTPBP2', 'H1FX', 'H1FX-AS1', 'HAUS4', 'HAUS7',
    'HBA1', 'HBA2', 'HBE1', 'HBG2', 'HES4', 'HIBADH', 'HIST1H1B',
    'HIST1H1E', 'HIST1H2AB', 'HIST1H2AH', 'HIST1H2AI', 'HIST1H2AJ',
    'HIST1H2BB', 'HIST1H2BJ', 'HIST1H2BK', 'HIST1H2BM', 'HIST1H2BO',
    'HIST1H3C', 'HIST1H3F', 'HIST1H3J', 'HIST1H4F', 'HIST1H4J',
    'HIST2H2AA3', 'HIST2H2AA4', 'HIST2H3D', 'HIST2H4B', 'HK3', 'HLA-C',
    'HLA-DQA2', 'HLA-DQB2', 'HLA-DRB5', 'HLA-G', 'HNRNPA1P48',
    'HOMER1', 'HOXB2', 'IDUA', 'IFITM1', 'IGFLR1', 'IGHG1', 'IGHG2',
    'IGHG3', 'IGHG4', 'IGHGP', 'IGKV3-20', 'IGKV4-1', 'IGLC2', 'IGLC3',
    'IGLC6', 'IGLC7', 'IGLV1-44', 'IGLV2-18', 'IGLV2-23', 'IGLV2-8',
    'IGSF10', 'IL10RB', 'IL12A', 'IL12A-AS1', 'ITIH5', 'KCNC1',
    'KCTD13', 'KMT2E-AS1', 'KRBA2', 'KRTCAP2', 'LAMC2', 'LARGE1',
    'LATS2-AS1', 'LDB1', 'LENG9', 'LGALS9', 'LGALS9C', 'LILRA5',
    'LINC00240', 'LINC00513', 'LINC01001', 'LINC01116', 'LINC01550',
    'LINC01641', 'LINC01659', 'LINC01800', 'LINC01871', 'LINC01943',
    'LINC02062', 'LINC02160', 'LINC02256', 'LINC02394', 'LINC02415',
    'LINC02422', 'LIPE', 'LIX1L-AS1', 'LPAR4', 'LPAR6', 'LRP2BP',
    'LRP5', 'LRRC26', 'LRRC7', 'LY6E', 'LY6G6F', 'LY6G6F-LY6G6D',
    'LYG2', 'MAD2L1BP', 'MAFIP', 'MAS1', 'MCEE', 'MDP1', 'MEG8',
    'METRN', 'MFSD14A', 'MIA-RAB4B', 'MIF', 'MIR133A1HG', 'MPRIP-AS1',
    'MPV17L', 'MRPL30', 'MRPL38', 'MRPL39', 'MRPL53', 'MRPS17',
    'MT-ATP8', 'MT-CO2', 'MT-ND4L', 'MTFP1', 'MTRNR2L1', 'MTRNR2L10',
    'MTRNR2L11', 'MTRNR2L12', 'MTRNR2L7', 'MTRNR2L8', 'MYCBP', 'MYH7B',
    'NACA2', 'NAV2-AS2', 'NBEAL1', 'NBPF12', 'NBPF14', 'NBPF20',
    'NCSTN', 'NDUFA13', 'NDUFA7', 'NDUFB3', 'NME2', 'NOMO1', 'NOMO2',
    'NPIPA1', 'NPIPB4', 'NPW', 'NRROS', 'NUTM2A-AS1', 'OR11G2',
    'OR2C3', 'OR4D9', 'OSBPL10-AS1', 'PABPC3', 'PADI4', 'PAGR1',
    'PANX2', 'PART1', 'PBX2', 'PDCD4-AS1', 'PEBP4', 'PGF', 'PIK3R2',
    'PKN1', 'PLCG2', 'PLPPR3', 'POLD4', 'POLR2J2', 'POLR2J3',
    'PPIAL4G', 'PRG2', 'PRKCD', 'PRRG3', 'PSMA2', 'PSMB10', 'PSMD9',
    'PTGDR2', 'PTPMT1', 'PWP2', 'RAB24', 'RAB38', 'RABAC1', 'RABGGTB',
    'RASD1', 'RAVER1', 'RDH5', 'RENBP', 'RGCC', 'RGS14', 'RGSL1',
    'RNASEH2A', 'RNASEK-C17orf49', 'RNF139-AS1', 'RPL17', 'RPL23A',
    'RPL36A', 'RPL39L', 'RPS10', 'RPS17', 'RPS18', 'RPS28', 'S100A1',
    'SAP30L-AS1', 'SCAND1', 'SCGB1C1', 'SCGB1C2', 'SCT', 'SEMA3E',
    'SEMA6B', 'SERF1B', 'SGSM3', 'SH2B2', 'SH3D19', 'SHISA8', 'SHPK',
    'SIDT1-AS1', 'SLC35B2', 'SLX1A', 'SMIM37', 'SNHG25', 'SNX15',
    'SP2-AS1', 'SPATC1L', 'SPDYC', 'SPDYE5', 'SPESP1', 'ST8SIA6-AS1',
    'STAG3', 'SUSD6', 'SVOPL', 'TBC1D3D', 'TCF23', 'TCF25', 'TEKT4',
    'TGFB1', 'THOC3', 'TJP2', 'TLCD2', 'TMC7', 'TMEM191B', 'TMEM216',
    'TMEM219', 'TMEM238', 'TOMM5', 'TONSL', 'TPSAB1', 'TPSB2', 'TPSD1',
    'TPTEP2-CSNK1E', 'TRAPPC5', 'TRBVB', 'TRGC1', 'TRGC2', 'TRPT1',
    'TSC22D4', 'TSHR', 'TUBA8', 'TVP23C-CDRT4', 'TWF2', 'TXNDC5',
    'U2AF1', 'U2AF1L5', 'UBE2V1', 'UCKL1', 'UGDH-AS1', 'UNC50',
    'UPK3BL1', 'USF1', 'USP19', 'WASHC1', 'WASIR1', 'WASIR2', 'WDR74',
    'XCL2', 'YIF1A', 'YIF1B', 'Z98884.1', 'ZBTB20-AS3', 'ZBTB7B',
    'ZC3H11A', 'ZDHHC8', 'ZFP36', 'ZNF414', 'ZNF670-ZNF695', 'ZNF865',
    'ZP3',
]
ULT_ILL_PROBLEMATIC_GENES = [
    *GENES_WITH_MORE_THAN_4_ULT_ILL_FOLD,
    *GENES_WITH_HIGH_ULT_ILL_LOG_RATIO_STD_AND_AT_LEAST_SOMETIME_ULT_ILL_DIFF_EXP,
]

# BATCHY_IN_ANY_STATE_SET_BY_KRUSKAL_GENES = [
#     'AC007952.4', 'AC092171.3', 'AC103591.3', 'AC245014.3', 'ACADVL', 'ADSL', 'AES', 'AL360012.1', 'AP2B1', 'ARL6IP1', 'ATP2C1', 'ATP5MD', 'ATP5ME', 'ATP6V0A1', 'ATP6V1G1', 'C1orf56', 'C6orf48', 'C8orf59', 'CALR', 'CAPNS1', 'CBX3', 'CCDC85B', 'CCND2', 'CCNL1', 'CD81', 'CD82', 'CDC42', 'CDV3', 'CETN3', 'CHCHD2', 'CHST12', 'CLIC1', 'COMMD6', 'CORO1A', 'COX16', 'COX6C', 'COX7A2', 'COX7B', 'DDAH2', 'DDX1', 'DDX21', 'DDX24', 'DDX3X', 'DDX5', 'DPYSL2', 'DUSP1', 'DYNC1I2', 'EEF1G', 'EGFL7', 'EGLN2', 'EIF2S2', 'EIF4B', 'EIF4H', 'ELP2', 'ENSA', 'FAM133B', 'FKBP1A', 'FNBP4', 'FXR1', 'GNA15', 'GNAI2', 'GNG11', 'GSN', 'H1FX', 'HCST', 'HIST1H1C', 'HIST1H1E', 'HIST1H2AC', 'HIST1H2BC', 'HIST1H4C', 'HIST1H4E', 'HIST1H4J', 'HIST2H2AC', 'HIST2H2BF', 'HM13', 'HNRNPA1', 'HNRNPA1P48', 'HNRNPA3', 'HNRNPDL', 'HNRNPL', 'HNRNPM', 'HSPE1', 'ICAM3', 'IGHA1', 'IST1', 'JUN', 'JUND', 'KLF6', 'KRT10', 'LAPTM5', 'LDHA', 'LRRC75A', 'LY6E', 'LYL1', 'LYZ', 'MAP3K7CL', 'MDH2', 'MDM4', 'MEF2C', 'MYL12A', 'NAA38', 'NAMPT', 'NARS', 'NBEAL1', 'NDUFA1', 'NDUFA13', 'NDUFA4', 'NDUFB1', 'NKTR', 'NMT1', 'PABPC4', 'PABPN1', 'PCBP2', 'PF4', 'PIH1D1', 'PKN2', 'PNRC1', 'PPBP', 'PPP1R14B', 'PRPF4B', 'PSMD14', 'PTCH2', 'PTPN6', 'R3HDM4', 'RANBP6', 'RASSF1', 'RBM39', 'RBM4', 'RBMX', 'RBX1', 'RGS18', 'RHOG', 'RNF145', 'RPL22L1', 'RPL36A', 'RPS17', 'S100A8', 'S100A9', 'SARAF', 'SCAND1', 'SEC61G', 'SELENOH', 'SEM1', 'SETD2', 'SF3B1', 'SFPQ', 'SHOC2', 'SMCHD1', 'SMIM26', 'SNHG25', 'SNHG9', 'SNRNP70', 'SNRPB2', 'SPCS2', 'SRRM1', 'SRSF10', 'SYF2', 'TBCB', 'TFPI', 'TMEM107', 'TMSB4X', 'TRIM8', 'TUBA1A', 'TUBB1', 'UBC', 'UBE2Q1', 'UQCR11', 'WDR45B', 'WDR74', 'YBX1', 'YWHAH', 'Z93241.1', 'ZBTB16', 'ZFC3H1', 'ZNF439', 
# ]
# UPREGULATED_IN_ANY_STATE_SET_IN_DELAYED_GENES = [
#     'AC013394.1','AC025159.1','AC026979.2','AC044849.1','AC095055.1','AC103591.3','ADAR','ADIPOR1','ADNP2','AFF1','AFF4','AKAP8L','AKIRIN1','AL355472.1','AL360012.1','AMD1','ANKRD12','ANKRD28','ANKRD36C','AP000547.3','ARF1','ARID2','ARID4A','ARID4B','ARL8A','ARMCX3','ARPC5L','ARPP19','ATAD2','ATF3','ATF4','ATF7IP2','ATG14','ATP2C1','ATP6V0A1','BAALC','BAZ2A','BCL3','BCLAF1','BEX2','BHLHE40','BIRC2','BNIP2','BNIP3L','BORCS7','BSDC1','BTG1','BTG2','BTG3','BUD31','BZW1','C16orf72','C16orf87','C18orf25','C6orf48','C9orf78','CCDC59','CCND2','CCNL1','CD69','CD82','CDC14A','CDC42','CDC42SE2','CDK11A','CDKN1B','CDKN2D','CDV3','CEP95','CGGBP1','CHD1','CHD2','CHD9','CHST11','CIR1','CIRBP','CKAP2','CLDND1','CLEC2B','CLIC1','CLINT1','CNOT1','CNST','COG3','COPS2','COQ10B','CREBRF','CSDE1','CSGALNACT2','CSNK1D','CTDSPL2','CTSW','CWC25','CYCS','CYLD','DAZAP1','DAZAP2','DCP1A','DDAH2','DDIT3','DDIT4','DDX21','DDX24','DDX3X','DDX5','DDX6','DEK','DENND3','DENND5A','DLD','DMXL1','DNAJB6','DNAJC7','DNTTIP2','DRAP1','DUSP10','DYNC1H1','DYNLL2','EAPP','EIF1','EIF1B','EIF2AK1','EIF4B','EIF4H','EIF5A','ELF1','ELOVL5','EMD','EMP3','ENSA','ENY2','EPB41L4A-AS1','EPC1','EREG','ETF1','ETFDH','EWSR1','EZH2','FAM133B','FAM160B1','FAM53C','FAM91A1','FBRS','FBXL20','FBXO7','FMNL1','FNBP4','FNIP1','FUBP1','FYTTD1',
#     ###'GATA2', # 240328: should we really mark GATA2 as lateral and noisy? and all of the others here???
#     'GCC2','GLUL','GNA15','GNAS','GNB1','GOLGA8B','GOLPH3','H1F0','H1FX','H2AFJ','H2AFV','H3F3B','HDGFL3','HEBP2','HIF1A','HIPK1','HIST1H1C','HIST1H1D','HIST1H1E','HIST1H2AK','HIST1H4E','HIST2H2AC','HLA-C','HMGA2','HMGB3','HNRNPA0','HNRNPA3','HNRNPC','HNRNPD','HNRNPDL','HNRNPH3','HNRNPL','HNRNPM','HPS4','HSF2','IER2','IER5','IGHM','IP6K1','IRF2BP2','IRGQ','ITM2B','ITM2C','ITSN2','JARID2','JMJD1C','JMJD6','JOSD1','JUN','JUNB','JUND','KAT7','KCMF1','KDM6A','KLF13','KLF3','KLF6','KLF9','KLHL24','KMT2C','KMT2E','KPNA3','KRR1','LAPTM4A','LAPTM5','LCOR','LDHA','LEPROT','LEPROTL1','LINC02573','LMAN1','LMO2','LSM14A','LYN','MAFG','MAN2A2','MAP1LC3B','MAP3K8','MAPK1IP1L','MAPKAPK2','MAX','MCL1','MDM4','MED13L','MED30','MEX3C','MIF','MIR181A1HG','MIR222HG','MIR4435-2HG','MOB1B','MPC2','MRPL22','MSI2','MTHFD2L','MTPN','MTRF1L','MTURN','MXD1','MXD4','N4BP2L1','NCOA7','NDUFA4','NDUFB4','NDUFS5','NECAP1','NFKB2','NFKBIA','NFKBIZ','NKTR','NORAD','NRIP1','NSD3','NT5C3A','NUP98','OCIAD2','OPA1','OPTN','OSER1','PABPN1','PCBP1','PCIF1','PCNP','PDIA3','PDRG1','PDZD8','PER1','PHACTR4','PHIP','PHTF1','PIK3IP1','PIK3R1','PIM3','PJA2','PKN2','PLAC8','PLCG2','PLEKHO1','PLIN2','PMAIP1','PNN','PNRC1','PNRC2','PPDPF','PPM1A','PPP1CB','PPP1R15A','PPP1R16B','PPP1R2','PPP2CA','PPP2R2D','PPP4R2','PRKAR1A','PRKD3','PRPF40A','PRPF4B','PRRC2C','PSMB1','PSMD7','PSME4','PTAR1','PTBP2','PTBP3','PTEN','PTP4A1','PTP4A2','PTTG1IP','R3HDM4','RAB11FIP1','RAB21','RAB5A','RANBP2','RASSF1','RB1CC1','RBBP6','RBM3','RBM33','RBM38','RBM39','RBM4','RC3H2','REL','RELA','RELB','REST','RFLNB','RHOG','RICTOR','RIN3','RIOK3','RIPK2','RIT1','RLF','RNF10','RNF103','RNF11','RNF114','RNF138','RNF145','RNF181','RNF19A','RNF7','RNMT','RNPS1','ROCK1','RPL22L1','RRBP1','RSBN1','RSL24D1','RSRC2','RTF1','SAMD8','SAP18','SAT1','SBDS','SCAND1','SDCBP','SDE2','SEC61B','SEC61G','SELENOK','SEPT2','SERINC1','SERTAD2','SETD2','SETD7','SF1','SF3B1','SFPQ','SH3GLB1','SHOC2','SIAH2','SINHCAF','SKIL','SLC38A2','SLU7','SMARCE1','SMCHD1','SMIM3','SNHG15','SNHG7','SNHG9','SNRNP70','SNRPA1','SNW1','SPAG9','SPTBN1','SQSTM1','SREK1','SRRM1','SRRM2','SRSF10','SRSF11','SRSF5','SSFA2','SSH2','STAM','STAT3','STAU1','STK17B','STK4','STRAP','STX12','SUMO3','SUZ12','SYAP1','SYF2','SYNGR1','TACC1','TAF1D','TAF9','TAOK3','TAPBP','TAPT1','TAX1BP1','TCF25','TEC','TESMIN','TGOLN2','TLK2','TM9SF3','TMEM123','TMEM167B','TMEM41B','TOMM20','TOPORS','TOX4','TP53INP1','TPM4','TRA2B','TRAF4','TRIM8','TSPYL2','TWISTNB','UBALD2','UBB','UBC','UBE2A','UBE2B','UBE2D2','UBE2D3','UBE2H','UBE2J1','UBQLN1','UBXN2A','UBXN4','UHRF2','UNK','UQCRFS1','USP12','USP14','USP22','VAV3','VIM','VPS26A','VPS37B','WAC','WBP11','WDR45B','WDR48','WHAMM','YBX3','YME1L1','YPEL5','YTHDC1','YWHAH','YWHAZ','ZBTB10','ZBTB21','ZC3H15','ZC3H7A','ZFC3H1','ZFR','ZNF292','ZNF326','ZNF555','ZNF639','ZNF711','ZNF738',
# ]

TOP2A_MKI67_GENE_CLUSTER = [
    # 230227: the following were in a module containing TOP2A and MKI67 (ordered by sub-modules):
    *['ASPM', 'BIRC5', 'CDKN3', 'CCNB2'],
    *['MKI67', 'TOP2A', 'POLQ', 'TTK', 'CDK1', 'CDCA5', 'DIAPH3'],
    *[
        'RRM2', 'TK1', 
        'UBE2C', # a well-known one. https://febs.onlinelibrary.wiley.com/doi/full/10.1111/febs.15134: "Ubiquitin-conjugating enzyme 2C (UBE2C) is a core ubiquitin-conjugating enzyme in the ubiquitin–proteasome system that promotes cell cycle progression" and "UBE2C, a ubiquitin-conjugating enzyme, accepts ubiquitin from E1, transfers it to specific anaphase-promoting complex/cyclosome (APC/C) substrates and catalyses lys-11- and lys-48-specific polyubiquitination, finally contributing to degradation of the APC/C substrates via the ubiquitin–proteasome system and promoting mitotic exit and cell cycle progression [13-15]. Most cancer cells are characterized by rapid proliferation and cell cycle dysregulation. Recent accumulating evidence suggests that UBE2C is one of a variety of genes that regulate cell cycle progression".
        'NCAPH', 'MELK',
    ],
    *['HJURP', 'DLGAP5', 'CDC20', 'CEP55', 'HIST1H3G', 'SHCBP1', 'NEK2', 'HMMR', 'CDCA8', 'NCAPG'],
    *['DEPDC1', 'CKAP2L', 'CCNA2', 'KIF2C', 'KIF4A', 'CIT', 'SPC24'],
    *['CENPA', 'CDCA2', 'GTSE1', 'FOXM1', 'SPAG5', 'CENPE', 'KIF23', 'TROAP', 'NUF2', 'KIFC1'],
]

Y_LINKED_GENE_NAMES = [
    # not exhaustive. just took genes that i already annotated as sex diff expressed and google says they are Y linked. can't just take all genes on Y because of pseudo-autosomal regions.
    'RPS4Y1', # https://www.proteinatlas.org/ENSG00000129824-RPS4Y1: "Ribosomal protein S4 Y-linked 1"

    'TMSB4Y', # https://www.proteinatlas.org/ENSG00000154620-TMSB4Y: "Thymosin beta 4 Y-linked"
    'UTY', # https://www.proteinatlas.org/ENSG00000183878-UTY: "Ubiquitously transcribed tetratricopeptide repeat containing, Y-linked" and "Male-specific histone demethylase that catalyzes trimethylated 'Lys-27' (H3K27me3) demethylation in histone H3. Has relatively low lysine demethylase activity."
    'DDX3Y', # https://www.proteinatlas.org/ENSG00000067048-DDX3Y: "DEAD-box helicase 3 Y-linked" and "Probable ATP-dependent RNA helicase. During immune response, may enhance IFNB1 expression via IRF3/IRF7 pathway (By similarity)."
    'ZFY', # https://www.proteinatlas.org/ENSG00000067646-ZFY: "Zinc finger protein Y-linked"
    'EIF1AY', # https://www.proteinatlas.org/ENSG00000198692-EIF1AY: "Eukaryotic translation initiation factor 1A Y-linked" and "Seems to be required for maximal rate of protein biosynthesis. Enhances ribosome dissociation into subunits and stabilizes the binding of the initiator Met-tRNA(I) to 40 S ribosomal subunits (By similarity)." and "This gene is located on the non-recombining region of the Y chromosome. It encodes a protein related to eukaryotic translation initiation factor 1A (EIF1A), which may function in stabilizing the binding of the initiator Met-tRNA to 40S ribosomal subunits."
    'KDM5D', # https://www.proteinatlas.org/ENSG00000012817-KDM5D: "Lysine demethylase 5D" and "Histone demethylase that specifically demethylates 'Lys-4' of histone H3, thereby playing a central role in histone code". https://pubmed.ncbi.nlm.nih.gov/32081420/: "The KDM5C and KDM5D genes, which encode H3K4 histone demethylases, are a surviving ancestral gene pair located on the X and Y chromosomes, respectively."
]

MT_RNR2_LIKE_GENE_NAMES = [
    'MTRNR2L10',
    'MTRNR2L11',
    'MTRNR2L12',
    'MTRNR2L1',
    'MTRNR2L13',
    'MTRNR2L6',
    'MTRNR2L3',
    'MTRNR2L5',
    'MTRNR2L7',
    'MTRNR2L4',
    'MTRNR2L8',
]

LATERAL_GENES_THAT_MIGHT_BE_MISSING_FROM_C_AD = [
    # 'TBC1D3D',
    # 'ATXN7-1', 'MATR3-1', 'TBCE-1',

    # # missing from only_N257_bm_06_10_22_1 (excluded due to zero UMIs in all droplets, i guess)
    # 'IGHV5-51', 'IGKV6-21', 'IGLV5-48', 'IGHV3-69-1', 'IGHVIII-44', 'IGHV3-19', 'IGHVII-74-1', 'IGHVII-44-2', 'IGLV1-36', 'IGLV10-54', 'IGHGP', 'IGHV1OR15-2', 'IGKV1-27', 'IGHV4-80', 'IGKV3D-11', 'IGHV3-15', 'IGKV1OR2-3', 'Z98884.1', 'IGKV1D-17', 'IGHV1-58', 'AC109597.2', 'IGHV3-66', 'IGHV3OR16-9', 'IGHV7-27', 'IGHVII-22-1', 'IGHV7-40', 'IGHV3OR16-7', 'IGKV2-29', 'IGLVIVOR22-1', 'IGHVIII-2-1', 'IGHV2-70', 'IGLC4', 'IGHV1-67', 'IGHV3-52', 'IGHV3-21', 'IGHJ4', 'IGLV3-10', 'IGHVIV-44-1', 'IGLVI-56', 'IGHVII-26-2', 'IGHV1-69', 'IGHV2OR16-5', 'IGHV7-4-1', 'IGKV2-24', 'IGLV2-23', 'IGHV4-61', 'AC007014.1', 'IGHV3-13', 'IGLV3-17', 'IGLV1-40', 'IGHV3-29', 'IGKV1D-8', 'IGLV7-35', 'IGKV6D-21', 'IGKV3-7', 'IGHVII-30-21', 'IGLC5', 'IGHV3-30-2', 'IGHV1OR15-3', 'IGKV3D-15', 'IGLVV-58', 'IGHEP1', 'IGKV1D-37', 'IGLV5-37', 'IGHV5-78', 'IGHV4-4', 'IGKV1-17', 'IGHV3-53', 'IGKV1D-33', 'IGHV3-62', 'DNAJB8-AS1', 'IGHV1OR16-4', 'IGLV4-60', 'IGKV2D-30', 'IGLV9-49', 'IGKV2-26', 'IGKV1-9', 'IGKV3D-20', 'IGLV3-19', 'MTRNR2L13', 'IGKV1-39', 'DNAJB8', 'IGKV2-30', 'IGHV3-42', 'IGHV4-28', 'IGHV1-14', 'IGLVIV-59', 'IGHV3OR16-6', 'IGHV3-64D', 'IGHV1-3', 'IGHV4-31', 'IGHV3-72', 'IGHV3-64', 'IGHJ6', 'IGLV3-12', 'IGKV1D-43', 'IGHV3-22', 'IGLV7-46', 'IGLV3-9', 'IGHV2-70D', 'IGLV3-16', 'IGHV1-24', 'IGKV2OR22-3', 'IGHV3-65', 'IGLV3-25', 'IGLV5-45', 'IGLV2-18', 'IGHV3-33', 'IGHV1-18', 'IGLV2-5', 'IGHV3-11', 'IGLV7-43', 'IGHV3-71', 'AC010883.1', 'IGLV1-50', 'SEMA3E', 'IGHVII-28-1', 'IGHV2-26', 'IGKV2D-29', 'IGLC6', 'IGKV2OR2-1', 'IGKV1D-13', 'IGLV3-27', 'IGKV2-10', 'IGLV1-41', 'IGLV1-44', 'IGKV1-12', 'IGLV4-69', 'IGHE', 'IGHV3-41', 'IGHV6-1', 'IGHV1OR15-4', 'IGLV3-26', 'IGHV3-48', 'IGLV6-57', 'TMSB4Y', 'LINC01116', 'IGKV1-6', 'IGHV3-38', 'IGLV8OR8-1', 'IGHVII-30-1', 'IGHV4-39', 'IGLV2-11', 'IGHVII-60-1', 'HLA-DQB1-AS1', 'PRRG3', 'IGLVI-70', 'IGHV1-68', 'IGKV1D-16', 'IGHV1-12', 'IGHV3-74', 'IGHV1-69-2', 'IGHV3-43', 'IGHV3-57', 'IGHV1-17', 'CXCL5', 'IGLV3-1', 'DDX3Y', 'IGKV2D-28', 'RPS4Y2', 'IGHVII-33-1', 'IGHV5-10-1', 'IGHV1OR15-6', 'IGKV1D-39', 'GPM6A', 'IGKV1-8', 'IGHV1-46', 'AC211486.2', 'IGHJ2P', 'IGLJ3', 'IGKV1OR2-108', 'KDM5D', 'IGKV1-37', 'IGKV7-3', 'IGKV2-18', 'IGKV2D-24', 'IGLV3-21', 'IGHV3-32', 'IGLV8-61', 'MTRNR2L5',

    # # missing from only_demux_06_06_22_2_ultima (excluded due to zero UMIs in all droplets, i guess)
    # 'SCT', 'IGKV5-2', 'IGLV1-47', 'DNAJB5-DT', 'IGHV4-59', 'AP002449.1', 'IGLV4-3', 'AC004556.1', 'IGHV3-79', 'IGHV3-75', 'AL133523.1', 'IGLV2-14', 'IGHV1-69D',
    
    # # missing from only_N328 (excluded due to zero UMIs in all droplets, i guess)
    # 'RPS4Y1', 'IGLC5', 'IGHV1-18', 'IGHV3-13', 'IGKV1D-37', 'IGHV4-55', 'IGLV2-8', 'IGHV3OR16-9', 'IGKV3-7', 'IGKV1-16', 'IGKV2D-29', 'IGHV2-70D', 'IGLV4-60', 'IGKV1-17', 'IGLV4-3', 'IGKV2OR2-1', 'IGKV1D-39', 'IGKV2D-24', 'IGLV2-18', 'IGHV4-80', 'IGHV3-16', 'IGLV1-36', 'IGHV3OR15-7', 'IGHJ4', 'IGKV2D-28', 'IGLV2-5', 'IGKV1-27', 'IGKV1OR2-108', 'IGHV3-52', 'IGHV4-39', 'IGKV1D-43', 'IGHVII-30-1', 'IGLV3-6', 'EIF1AY', 'IGKV2-30', 'IGLVI-56', 'IGKV2D-30', 'IGHVII-28-1', 'IGLC4', 'IGHV4-28', 'IGHV3-19', 'IGLV3-9', 'IGKV1-8', 'IGHV4-31', 'PRRG3', 'IGHV3-64', 'IGLV2-23', 'IGHV3-73', 'IGLJ3', 'IGKV1D-16', 'IGLV3-10', 'IGHV3-15', 'IGHG1', 'IGKV3D-11', 'IGHV1-67', 'IGHV3-29', 'IGHV2-26', 'IGKV1OR22-1', 'DNAJB8', 'IGKV1-12', 'IGKV2D-26', 'IGHV3-22', 'IGLV7-46', 'IGHV5-51', 'IGKV3-11', 'IGHV3-25', 'IGLON5', 'IGLV7-43', 'IGLV1-50', 'IGHVII-65-1', 'IGLVI-70', 'IGHV3-21', 'RPS4Y2', 'IGLV8-61', 'IGLVIVOR22-1', 'IGKV2-29', 'IGHV2-70', 'IGLV1-40', 'IGKV2D-18', 'IGLV10-54', 'LINC02160', 'IGKV1D-13', 'IGKV1D-8', 'IGHVII-60-1', 'IGHVIII-2-1', 'IGHV1OR15-6', 'IGHV3-43', 'IGLV6-57', 'IGKV1D-17', 'IGLC7', 'IGHV3-11', 'IGHV6-1', 'IGLV3-25', 'IGKV3D-15', 'IGHG2', 'IGHEP1', 'IGLV1-44', 'IGHV1-45', 'IGHV3-48', 'IGKV1-9', 'IGKV1-37', 'IGLV3-21', 'IGLV5-37', 'AC007014.1', 'IGHV3-41', 'IGHV4-4', 'IGKV1-6', 'IGLV4-69', 'IGLV5-45', 'IGHV3-30', 'IGHVII-43-1', 'TMEM191B', 'IGHVII-22-1', 'IGHV4-61', 'IGKV1D-33', 'DNAJB8-AS1', 'IGKV1-5', 'IGLL5', 'IGHV7-4-1', 'IGHV1-12', 'IGHV3-64D', 'IGKV2-26', 'IGHV5-10-1', 'IGLV2-11', 'IGHV3-53', 'IGHV3-71', 'IGKV5-2', 'IGHV3-23', 'IGHV4-59', 'IGLV1-41', 'IGKV1D-12', 'IGKV3-15', 'IGHV1-46', 'IGLV3-27', 'IGLV3-22', 'SEMA3E', 'IGLV5-48', 'IGHV1-69D', 'IGKV1-39', 'IGHV3OR16-8', 'IGKV6-21', 'AC211486.2', 'IGHGP', 'IGKV2OR22-3', 'IGLV3-16', 'IGHV7-27', 'Z98884.1', 'KDM5D', 'IGHV1-24', 'IGLV1-47', 'IGLV2-14', 'IGLV9-49', 'IGLC6', 'IGHVIII-44', 'MTRNR2L13', 'IGHV3-32', 'IGHV3OR16-7', 'IGKV3D-20', 'IGLV3-19', 'IGHE', 'IGKV2-28',

    
    *MT_RNR2_LIKE_GENE_NAMES, 'NEAT1', 'XIST',
    
    # missing from '/dummy/dummy/dummy/raid/mds/230706_healthy_ill_with_24_36_48/intermediate_output/raw_without_excluded_genes_and_cells.h5ad'
    'MTRNR2L5',
    'IGHJ6', 'IGHV3-30-2', 'IGHVII-30-21', 'IGKV7-3', 'IGLV10-67', 'IGLV3-26', 'IGHV1OR16-1', 'IGLV2-34', 'IGLV3-2', 'IGHV3-6', 'IGHV7-40', 'IGKV2-4', 'IGLVIV-59', 'IGHV7-81', 'IGLV3-12', 'IGKV3-34', 'IGHV3-57', 'IGKV2-18', 'IGKV2-10', 'IGKV6D-21', 'IGHV3-42', 'IGHV2OR16-5', 'IGKV6D-41', 'IGHV1-68', 'IGLV3-17', 'IGLV7-35', 'IGHV1OR21-1', 'IGHV1-17', 'IGKV1OR2-3', 'IGKV1-13', 'IGHJ2P', 'IGHV3OR16-10', 'IGLV11-55', 'IGLV8OR8-1', 'IGHVIII-82', 'IGKV2OR2-2', 'IGKV1OR1-1', 'IGHV3-60', 'IGHVII-33-1', 'IGHV1OR16-4', 'IGKV2-23', 'IGKV1D-35', 'IGKV1D-27', 'IGHVII-44-2', 'IGKV1OR-3', 'IGHV1OR15-4', 'IGHV3-37', 'IGHV1-14', 'IGHV3-63', 'IGHV3OR16-6', 'IGHVIV-44-1', 'IGHV1OR15-2', 'IGHV1OR15-3', 'IGHVII-74-1', 'IGHVII-26-2', 'IGKV1-35', 'IGHV3-62', 'IGHV3OR16-13', 'IGHV3-33-2', 'IGHVII-78-1', 'IGLVV-58',

    # were missing from '/dummy/dummy/dummy/raid/mds/all_mpn_exps/intermediate_output/raw_without_excluded_genes_and_cells.h5ad' on 240123:
    'TRAV8-7', 'TRAV11', 'IGHV3OR16-8', 'TRBV12-5', 'IGHV3-25', 'IGHJ4', 'TRBV5-3', 'TRGV11', 'DNAJB8', 'TRBV7-1', 'TRBV20OR9-2', 'TRBV11-1',
    
    # were missing from '/dummy/dummy/dummy/raid/mds/new_N280_exps/intermediate_output/raw_without_excluded_genes_and_cells.h5ad' on 240129:
    'IGHV4-80', 'TRBV21-1', 'TRBV10-3', 'TRBV6-7', 'TRBV5-2', 'TRBV21OR9-2', 'TRBV12-2', 'TRAV10', 'TRBV29OR9-2', 'TRBV12-1', 'IGHV3-71', 'TRAV40', 'IGHVII-22-1', 'TRBV8-2', 'IGHV1OR15-6', 'TRAV12-1', 'DNAJB8-AS1', 'IGHV3-35', 'IGHVII-28-1', 'TRAV38-1', 'IGHVII-65-1', 'IGHVII-43-1', 'IGKV1D-17', 'IGHV1-67', 'TRAV18', 'TRAV8-1', 'IGHV3-52', 'IGKV1OR22-1', 'TRBV26', 'IGHVIII-2-1', 'TRBV6-1', 'TRBV1', 'TRBVA', 'IGKV1OR2-108', 'TRBV13', 'TRBV19', 'IGHV3-16', 'IGHV5-10-1', 'TRBV5-7', 'TRBV22-1', 'TRBV7-6', 'TRAV24', 'TRBV2', 'TRBV17', 'TRBV24-1', 'TRBV10-1', 'TRBV5-4', 'IGLV3-22', 'TRBV23OR9-2', 'TRBV25-1', 'TRBV11-2', 'IGLVIVOR22-1', 'IGHVII-30-1', 'TRGVB', 'IGHV3OR16-7', 'TRBV5-5', 'IGKV2-28', 'TRAV31', 'IGKV1D-37', 'TRBV6-4', 'IGHV3-29', 'IGHV3-65', 'TRAV2', 'TRAV28', 'TRAV30', 'TRBV7-7', 'IGKV1D-33', 'IGHV3-32', 'TRAV21', 'IGHV3-19', 'IGLVI-56', 'IGHV3OR16-9', 'TRAV9-2', 'TRBV8-1', 'TRBV4-2', 'TRBV11-3', 'TRAV37', 'IGHVIII-44', 'TRBV6-8', 'IGHV3-76', 'TRAV32', 'IGHV3-79', 'TRGV1', 'TRGVA', 'IGKV2OR2-1', 'TRAV1-1', 'IGHV1-12', 'TRAV34', 'TRAV41', 'IGHV3-64D', 'IGHV3-41', 'TRAV33', 'TRAV38-2DV8', 'TRBV10-2', 'IGKV2D-18', 'IGHV7-27', 'IGHV3-75',

    # were missing from '/dummy/dummy/dummy/raid/mds/240129_exps_with_donors_with_any_sample_that_waited/intermediate_output/raw_without_excluded_genes_and_cells.h5ad' on 240129:
    'TRAV26-2',

    # were missing from '/dummy/dummy/dummy/raid/mds/240222_pb_exps_potentially_containing_mpn_ult/intermediate_output/raw_without_excluded_genes_and_cells.h5ad' on 240228:
    'IGLV5-48', 'TRDJ1', 'IGLV3-6', 'IGHEP1', 'IGKV1-37', 'IGHV3OR15-7', 'IGKV1D-12', 'IGLV5-37',
    # were missing from '/dummy/dummy/dummy/raid/mds/240222_pb_exps_potentially_containing_mpn_ult/intermediate_output/raw_without_excluded_genes_and_cells.h5ad' on 240311:
    'TRAV16', 'IGLC4', 'IGHV4-61', 'TRAV8-4', 'IGLV4-3', 'TRAV1-2', 'TRAV27', 'TRAV12-2', 'TRBV30',
    # were missing from '/dummy/dummy/dummy/raid/mds/240222_pb_exps_potentially_containing_mpn_illu/intermediate_output/raw_without_excluded_genes_and_cells.h5ad' on 240228:
    'IGKV3D-15', 'TRBV15',
]

LATERAL_GENE_DICT = {
        'mt_rnr2_like': MT_RNR2_LIKE_GENE_NAMES,

        'misc': [
            # TODO: cleanup these misc lateral genes. maybe some of them should not be here.


            # TODO: mark AHNAK as lateral?
            # 221210: on the first metacell run, AHNAK was among the genes with highest 'feature_gene', and it looks like it is expressed in most of the metacells. in the umap it seems like it is enriched and many different places, and considering also its presumable function, which sounds relatively lateral, decided to mark it as lateral.
            
            # # AHNAK and genes correlated with AHNAK, but badly.
            # 'AHNAK', # https://www.proteinatlas.org/ENSG00000124942-AHNAK: "nucleoprotein" and "Tissue profilei	Ubiquitous cytoplasmic expression." and "Tissue specificityi	Low tissue specificity" and "Immune cell specificityi	Low immune cell specificity" and "structural scaffold protein" and "The encoded protein may play a role in such diverse processes as blood-brain barrier formation, cell structure and migration, cardiac calcium channel regulation, and tumor metastasis."
            # 'TSPAN2', 
            # 'ANXA1', # very strong. correlated with sum of all non-strong ones, but the MC scatter is a bit ugly. 
            # 'TAGLN2', # strong. correlated with sum of all non-strong ones, but the MC scatter is a bit ugly


            # TODO: should we do something about the distribution of TPT1?
            'TPT1', # looks a bit batchy - /dummy/dummy/dummy/tanay_group/mds/figs/221227_downsampled_TPT1_per_exp.png. more importantly, high relative variance, high expression, and involved in apoptosis and cell division. # https://www.proteinatlas.org/ENSG00000133112-TPT1: "This gene encodes a protein that is a regulator of cellular growth and proliferation. [...] The encoded protein is involved in a variety of cellular pathways, including apoptosis, protein synthesis and cell division. It binds to and stabilizes microtubules, and removal of this protein through phosphorylation is required for progression through mitotic and meiotic cell divisions."
            # # https://www.proteinatlas.org/ENSG00000133112-TPT1: "This gene encodes a protein that is a regulator of cellular growth and proliferation. Its mRNA is highly structured and contains an oligopyrimidine tract (5'-TOP) in its 5' untranslated region that functions to repress its translation under quiescent conditions. The encoded protein is involved in a variety of cellular pathways, including apoptosis, protein synthesis and cell division. It binds to and stabilizes microtubules, and removal of this protein through phosphorylation is required for progression through mitotic and meiotic cell divisions."
            # 'TPT1', # mc.pl.relate_genes() clustered it in a strong module of ribosomal protein genes (it was the only non ribosomal protein gene in the module). (the module was: ['RPS8', 'RPL32', 'RPL34', 'RPS3A', 'RPS23', 'RPS14', 'RPS18', 'RPS12', 'RPS4X', 'RPL39', 'RPL10', 'RPS6', 'RPL12', 'RPS3', 'RPS24', 'TPT1', 'RPLP1', 'RPS2', 'RPL13', 'RPL26', 'RPS19', 'RPL3'])

            # 221210: TUBA1A looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
            'TUBA1A', # https://www.proteinatlas.org/ENSG00000167552-TUBA1A: "Tubulin alpha 1a" and "Tubulin is the major constituent of microtubules. It binds two moles of GTP, one at an exchangeable site on the beta chain and one at a non-exchangeable site on the alpha chain."

            'TUBA1B', # https://www.proteinatlas.org/ENSG00000123416-TUBA1B: "Tubulin alpha 1b" and "Tubulin is the major constituent of microtubules."


            # 'POLR2J3', # https://www.proteinatlas.org/ENSG00000168255-POLR2J3: "RNA polymerase II subunit J3"

            'HMGB1', # https://www.proteinatlas.org/ENSG00000189403-HMGB1: "High mobility group box 1" and "Multifunctional redox sensitive protein with various roles in different cellular compartments. In the nucleus is one of the major chromatin-associated non-histone proteins and acts as a DNA chaperone involved in replication, transcription, chromatin remodeling, V(D)J recombination, DNA repair and genome stability 1. Proposed to be an universal biosensor for nucleic acids. Promotes host inflammatory response to sterile and infectious signals and is involved in the coordination and integration of innate and adaptive immune responses. In the cytoplasm functions as sensor and/or chaperone for immunogenic nucleic acids implicating the activation of TLR9-mediated immune responses, and mediates autophagy. Acts as danger associated molecular pattern (DAMP) molecule that amplifies immune responses during tissue injury"

            # # 221210: CALCOCO2 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
            # 'CALCOCO2', # https://www.proteinatlas.org/ENSG00000136436-CALCOCO2: "Calcium binding and coiled-coil domain 2" and "Acts as an effector protein of galectin-sensed membrane damage that restricts the proliferation of infecting pathogens such as Salmonella typhimurium upon entry into the cytosol by targeting LGALS8-associated bacteria for autophagy"

            # # 221210: SAMD9 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
            # 'SAMD9', # https://www.proteinatlas.org/ENSG00000205413-SAMD9: "Sterile alpha motif domain containing 9" and "May play a role in the inflammatory response to tissue injury and the control of extra-osseous calcification, acting as a downstream target of TNF-alpha signaling." and "localizes to the cytoplasm and may play a role in regulating cell proliferation and apoptosis."

            # # decided to not mark SAMD9L as lateral at this point because https://www.proteinatlas.org/ENSG00000177409-SAMD9L says "Naturally occurring mutations in this gene are associated with myeloid disorders such as juvenile myelomonocytic leukemia, acute myeloid leukemia, and myelodysplastic syndrome."
            # # 'SAMD9L',

            # # 221210: HLA-DPB1 and HLA-DRB1 don't look very nice in mcview, especially compared to how they look in nimrodra's mcview. TODO: mark theme as lateral?
            # # similarly CD74, which is in a highly correlated with multiple HLA genes.
            # # similarly CD52?? it is highly correlated with multiple HLA genes, but its role seems to be unknown.
            # # 'HLA-DPB1', # https://www.proteinatlas.org/ENSG00000223865-HLA-DPB1: "Major histocompatibility complex, class II, DP beta 1"
            # # 'HLA-DRB1', # https://www.proteinatlas.org/ENSG00000196126-HLA-DRB1: "Major histocompatibility complex, class II, DR beta 1"
            # # 'CD52', #
            # # 'HLA-DQB1', # https://www.proteinatlas.org/ENSG00000179344-HLA-DQB1
            # # 'HLA-DRA', # https://www.proteinatlas.org/ENSG00000204287-HLA-DRA
            # # 'HLA-DPA1', # https://www.proteinatlas.org/ENSG00000231389-HLA-DPA1
            # # 'CD74', # 221210: even though CD74 is highly correlated with multiple HLA genes, it looks quite nice in my mcview, so not marking it as lateral currently.

            # # 'HLA-C',
            # # 'HLA-DMB',
            # # 'HLA-A',
            # # 'HLA-E',
            # # 'HLA-B',
            # # 'HLA-DQA1',

            # # 221210: HLA-DRB5 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
            # # 'HLA-DRB5', # https://www.proteinatlas.org/ENSG00000198502-HLA-DRB5

            # # 221210: LINC01857 looked quite salt and pepper in the enrichment view in mcview, but as its function is unclear, I decided to maybe mark it as lateral later
            # # 'LINC01857', # perhaps relevant: https://www.nature.com/articles/s41417-020-00267-4: "LncRNA LINC01857 promotes cell growth and diminishes apoptosis via PI3K/mTOR pathway and EMT process by regulating miR-141-3p/MAP4K4 axis in diffuse large B-cell lymphoma"



            # # 221210: IGHMBP2 looked very salt and pepper in the enrichment view in mcview, and due to its presumable function (helicase - sounds not specific, in addition to data supporting it is not specific), decided to mark it as lateral.
            # 'IGHMBP2', # https://www.proteinatlas.org/ENSG00000132740-IGHMBP2: "5' to 3' helicase that unwinds RNA and DNA duplexes in an ATP-dependent reaction 1, 2, 3. Specific to 5'-phosphorylated single-stranded guanine-rich sequences 4, 5. May play a role in RNA metabolism, ribosome biogenesis or initiation of translation 6, 7. May play a role in regulation of transcription (By similarity). Interacts with tRNA-Tyr"

            # 221210: SMCHD1 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
            'SMCHD1', # https://www.proteinatlas.org/ENSG00000101596-SMCHD1: "Structural maintenance of chromosomes flexible hinge domain containing 1" and "Non-canonical member of the structural maintenance of chromosomes (SMC) protein family that plays a key role in epigenetic silencing by regulating chromatin architecture (By similarity). Promotes heterochromatin formation in both autosomes and chromosome X, probably by mediating the merge of chromatin compartments (By similarity). Plays a key role in chromosome X inactivation in females by promoting the spreading of heterochromatin 1. Recruited to inactivated chromosome X by Xist RNA and acts by mediating the merge of chromatin compartments: promotes random chromatin interactions that span the boundaries of existing structures, leading to create a compartment-less architecture typical of inactivated chromosome X (By similarity). Required to facilitate Xist RNA spreading (By similarity)."

            



            # 221210: HNRNPA1 seems almost ubiquitously expressed, and maybe a little bit salt and pepper. it is correlated with HSP90AB1, and its function sounds like it should be relatively lateral. not sure. also, highly expressed. also, looked a bit batchy.
            
            
            # 'HNRNPA1', # https://www.proteinatlas.org/ENSG00000135486-HNRNPA1: "Heterogeneous nuclear ribonucleoprotein A1" and "This gene encodes a member of a family of ubiquitously expressed heterogeneous nuclear ribonucleoproteins (hnRNPs), which are RNA-binding proteins that associate with pre-mRNAs in the nucleus and influence pre-mRNA processing, as well as other aspects of mRNA metabolism and transport. The protein encoded by this gene is one of the most abundant core proteins of hnRNP complexes and plays a key role in the regulation of alternative splicing."

            
            # 'HNRNPH1', # https://www.proteinatlas.org/ENSG00000169045-HNRNPH1: "Heterogeneous nuclear ribonucleoprotein H1" and "This gene encodes a member of a subfamily of ubiquitously expressed heterogeneous nuclear ribonucleoproteins (hnRNPs). The hnRNPs are RNA binding proteins that complex with heterogeneous nuclear RNA. These proteins are associated with pre-mRNAs in the nucleus and appear to influence pre-mRNA processing and other aspects of mRNA metabolism and transport."

            
            'ACTG1', # https://www.proteinatlas.org/ENSG00000184009-ACTG1: "Actin gamma 1" and "The beta and gamma actins co-exist in most cell types as components of the cytoskeleton and as mediators of internal cell motility. Actin gamma 1, encoded by this gene, is a cytoplasmic actin found in all cell types."

            # 221210: ACTB is a highly expressed gene and a strong feature gene currently (IIUC). also, its presumable function seems relatively lateral. also, looked a bit batchy.
            'ACTB', # https://www.proteinatlas.org/ENSG00000075624-ACTB: "Actin beta" and "The encoded protein is a major constituent of the contractile apparatus and one of the two nonmuscle cytoskeletal actins that are ubiquitously expressed."

            
            # 'PTMA', # https://www.proteinatlas.org/ENSG00000187514-PTMA: "Prothymosin alpha" and "Prothymosin alpha may mediate immune function by conferring resistance to certain opportunistic infections"

            # 221210: EEF1A1 is an extremely highly expressed gene and a strong feature gene currently (IIUC), and its presumable function sounds relatively lateral. also, looks a bit batchy.
            'EEF1A1', # https://www.proteinatlas.org/ENSG00000156508-EEF1A1: "Eukaryotic translation elongation factor 1 alpha 1" and "This protein promotes the GTP-dependent binding of aminoacyl-tRNA to the A-site of ribosomes during protein biosynthesis. Plays a role in the positive regulation of IFNG transcription in T-helper 1 cells as part of an IFNG promoter-binding complex with TXK and PARP1"
            

            # # 221210: GAPDH is a highly expressed gene and a strong feature gene currently (IIUC). also, its presumable function seems relatively lateral. doesn't look very batchy, so i think there isn't enough to mark it as lateral.
            # 'GAPDH', # https://www.proteinatlas.org/ENSG00000111640-GAPDH: "Glyceraldehyde-3-phosphate dehydrogenase" and "Has both glyceraldehyde-3-phosphate dehydrogenase and nitrosylase activities, thereby playing a role in glycolysis and nuclear functions, respectively 1, 2. Glyceraldehyde-3-phosphate dehydrogenase is a key enzyme in glycolysis that catalyzes the first step of the pathway by converting D-glyceraldehyde 3-phosphate (G3P) into 3-phospho-D-glyceroyl phosphate 3, 4. Modulates the organization and assembly of the cytoskeleton (By similarity). Facilitates the CHP1-dependent microtubule and membrane associations through its ability to stimulate the binding of CHP1 to microtubules (By similarity)."

            # 221210: TMSB10 is a highly expressed gene and a strong feature gene currently (IIUC). also, its presumable function seems relatively lateral. and it is in https://github.com/tanaylab/metacells/blob/master/vignettes/Metacells_Vignette.ipynb as one of the well known genes to exclude. and maybe looks a bit batchy. oh well.
            'TMSB10', # https://www.proteinatlas.org/ENSG00000034510-TMSB10: "Thymosin beta 10" and "Plays an important role in the organization of the cytoskeleton. Binds to and sequesters actin monomers (G actin) and therefore inhibits actin polymerization (By similarity)"

            
            # 'HNRNPA2B1', # https://www.proteinatlas.org/ENSG00000122566-HNRNPA2B1: "Heterogeneous nuclear ribonucleoprotein A2/B1" and "Heterogeneous nuclear ribonucleoprotein (hnRNP) that associates with nascent pre-mRNAs, packaging them into hnRNP particles. The hnRNP particle arrangement on nascent hnRNA is non-random and sequence-dependent and serves to condense and stabilize the transcripts and minimize tangling and knotting. Packaging plays a role in various processes such as transcription, pre-mRNA processing, RNA nuclear export, subcellular location, mRNA translation and stability of mature mRNAs"

            'KCNQ1OT1', 

            *GENES_EXCLUDED_IN_THE_METACELL_VIGNETTE,

            'IGLV2-14', # extremely correlated with IGLL5. https://www.proteinatlas.org/ENSG00000211666-IGLV2-14: "Immunoglobulin lambda variable 2-14" and "V region of the variable domain of immunoglobulin light chains that participates in the antigen recognition"



            # TODO: 'CTSS', 'MPEG1', 'FCN1', 'NCF2', 'MNDA', 'FGL2', # at least one of these has something suspicious going on in demux_09_02_22_1??

            # 'SAMHD1', # differentially expressed in GMP-L between demux_03_03_22_1 and demux_03_03_22_2, but the pvalue isnt crazy. maybe should be lateral anyway, but i am not sure. https://www.proteinatlas.org/ENSG00000101347-SAMHD1: "SAM and HD domain containing deoxynucleoside triphosphate triphosphohydrolase 1" and "Protein that acts both as a host restriction factor involved in defense response to virus and as a regulator of DNA end resection at stalled replication forks 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12. Has deoxynucleoside triphosphate (dNTPase) activity, which is required to restrict infection by viruses, such as HIV-1: dNTPase activity reduces cellular dNTP levels to levels too low for retroviral reverse transcription to occur, blocking early-stage virus replication in dendritic and other myeloid cells 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25. Likewise, suppresses LINE-1 retrotransposon activity"

            'MACF1', 

            'BPTF', 

            'PNISR', 

            'KMT2A', 

            'PRRC2C', 

            'HNRNPU', 

            'EIF2AK1', 

            'ZNF385D', # looked batchy in the kolmogorov smirnov tests when ignoring one donor of the experiment. https://www.proteinatlas.org/ENSG00000151789-ZNF385D: "Zinc finger protein 385D"
            
            'TUBB1', 

        ],

        'mitochondria_related': [
            # 230118: NDUFA13 looked very salt and pepper.
            'NDUFA13', # https://www.proteinatlas.org/ENSG00000186010-NDUFA13: "NADH:ubiquinone oxidoreductase subunit A13" and "Accessory subunit of the mitochondrial membrane respiratory chain NADH dehydrogenase (Complex I), that is believed not to be involved in catalysis 1. Complex I functions in the transfer of electrons from NADH to the respiratory chain. The immediate electron acceptor for the enzyme is believed to be ubiquinone 2. Involved in the interferon/all-trans-retinoic acid (IFN/RA) induced cell death. This apoptotic activity is inhibited by interaction with viral IRF1."

            # 230118: NDUFA6 looked salt and pepper.
            'NDUFA6', # https://www.proteinatlas.org/ENSG00000184983-NDUFA6: "NADH:ubiquinone oxidoreductase subunit A6" and "Accessory subunit of the mitochondrial membrane respiratory chain NADH dehydrogenase (Complex I), that is believed to be not involved in catalysis. Required for proper complex I assembly 1. Complex I functions in the transfer of electrons from NADH to the respiratory chain."

            'ATP5F1E', # looked batchy in the kolmogorov smirnov tests when ignoring one donor of the experiment. also, its function seems relatively lateral. https://www.proteinatlas.org/ENSG00000124172-ATP5F1E: "ATP synthase F1 subunit epsilon" and "Mitochondrial membrane ATP synthase (F(1)F(0) ATP synthase or Complex V) produces ATP from ADP in the presence of a proton gradient across the membrane which is generated by electron transport complexes of the respiratory chain."

            'COX7C', 

        ],

        'uncharacterized': [
            # TODO: is it actually ok to mark uncharacterized genes as lateral??? i guess the excuse is that we don't expect to identify a novel important gene here..

            
            

            # 'AL034397.3', # not marking as lateral because it is somewhat characterized, it seems. https://www.genecards.org/cgi-bin/carddisp.pl?gene=MIR223HG: "MIR223HG (MIR223 Host Gene) is an RNA Gene, and is affiliated with the lncRNA class. Diseases associated with MIR223HG include Leukemia, Acute Myeloid. Among its related pathways are miRNAs involvement in the immune response in sepsis."

            
        ],
    }
NOISY_GENE_DICT = {
    'sex_diff_expressed': [
        'XIST', # 221210: salt and pepper on the first metacell run, which could have been predicted by its differential expression between females and males.
        *Y_LINKED_GENE_NAMES,
        'RPS4X',
        'JPX', # 240403: seems higher in females, but not all of them (as in XIST). https://www.genecards.org/cgi-bin/carddisp.pl?gene=JPX: "JPX is a nonprotein-coding RNA transcribed from a gene within the X-inactivation center (XIC; MIM 314670) that appears to participate in X chromosome inactivation"
    ],
    'cell_cycle_related': [
        # DNA replication and cell division
        # 221210: each of MKI67, TOP2A, PCNA, ASPM, H2AFZ, MCM10 looked very/quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        # 230828: moved this into noisy genes because in the current healthy BM model, TOP2A, MKI67, etc, are among the genes with the highest number of metacells in which they were found to have significant inner fold.
        'MKI67',
        'TOP2A',
        *TOP2A_MKI67_GENE_CLUSTER,

        'PCNA',
        'POLD4', # https://www.proteinatlas.org/ENSG00000175482-POLD4: "DNA polymerase delta 4, accessory subunit" and "As a component of the tetrameric DNA polymerase delta complex (Pol-delta4), plays a role in high fidelity genome replication and repair. Within this complex, increases the rate of DNA synthesis and decreases fidelity by regulating POLD1 polymerase and proofreading 3' to 5' exonuclease activity 1, 2, 3. Pol-delta4 participates in Okazaki fragment processing, through both the short flap pathway, as well as a nick translation system 4. Under conditions of DNA replication stress, required for the repair of broken replication forks through break-induced replication (BIR), a mechanism that may induce segmental genomic duplications of up to 200 kb"
        'ASPM', # https://www.proteinatlas.org/ENSG00000066279-ASPM: "Assembly factor for spindle microtubules" and "Involved in mitotic spindle regulation and coordination of mitotic processes." and "This gene is the human ortholog of the Drosophila melanogaster 'abnormal spindle' gene (asp), which is essential for normal mitotic spindle function in embryonic neuroblasts. Studies in mouse also suggest a role of this gene in mitotic spindle regulation, with a preferential role in regulating neurogenesis."
        'H2AFZ', # https://www.proteinatlas.org/ENSG00000164032-H2AZ1: "H2A.Z variant histone 1" and "Gene namei H2AZ1 (H2A.Z, H2AFZ, H2AZ)"
        'MCM10', # https://www.proteinatlas.org/ENSG00000065328-MCM10: "Minichromosome maintenance 10 replication initiation factor" and "Acts as a replication initiation factor that brings together the MCM2-7 helicase and the DNA polymerase alpha/primase complex in order to initiate DNA replication. Additionally, plays a role in preventing DNA damage during replication."
        'CENPE', # highly correlated with TOP2A. https://www.proteinatlas.org/ENSG00000138778-CENPE: "Centromere protein E" and "Microtubule plus-end-directed kinetochore motor which plays an important role in chromosome congression, microtubule-kinetochore conjugation and spindle assembly checkpoint activation."
        'CENPF', # highly correlated with TOP2A. https://www.proteinatlas.org/ENSG00000117724-CENPF: "Centromere protein F" and "Required for kinetochore function and chromosome segregation in mitosis."
        'CENPU', # https://www.proteinatlas.org/ENSG00000151725-CENPU: "Centromere protein U" and "Component of the CENPA-NAC (nucleosome-associated) complex, a complex that plays a central role in assembly of kinetochore proteins, mitotic progression and chromosome segregation."
        'CENPA', # highly correlated with TOP2A. https://www.proteinatlas.org/ENSG00000115163-CENPA: "Centromere protein A" and "Required for recruitment and assembly of kinetochore proteins, and as a consequence required for progress through mitosis, chromosome segregation and cytokinesis", though also "The protein is a replication-independent histone that is a member of the histone H3 family."
        'CENPM', # https://www.proteinatlas.org/ENSG00000100162-CENPM: "Centromere protein M" and "Component of the CENPA-NAC (nucleosome-associated) complex, a complex that plays a central role in assembly of kinetochore proteins, mitotic progression and chromosome segregation"
        'CENPJ', # https://www.proteinatlas.org/ENSG00000151849-CENPJ: "Centromere protein J" and "Plays an important role in cell division and centrosome function by participating in centriole duplication"
        'CENPH', # https://www.proteinatlas.org/ENSG00000153044-CENPH: "Centromere protein H" and "Component of the CENPA-NAC (nucleosome-associated) complex, a complex that plays a central role in assembly of kinetochore proteins, mitotic progression and chromosome segregation".
        'SEPT7', # https://www.proteinatlas.org/ENSG00000122545-SEPTIN7: "Septin 7" and "Filament-forming cytoskeletal GTPase. Required for normal organization of the actin cytoskeleton. Required for normal progress through mitosis. Involved in cytokinesis. Required for normal association of CENPE with the kinetochore. Plays a role in ciliogenesis and collective cell movements".
        'KIF2C', # highly correlated with TOP2A. https://www.proteinatlas.org/ENSG00000142945-KIF2C: "Kinesin family member 2C" and "In complex with KIF18B, constitutes the major microtubule plus-end depolymerizing activity in mitotic cells 1. Regulates the turnover of microtubules at the kinetochore and functions in chromosome segregation during mitosis"
        'HMGB2', # https://www.proteinatlas.org/ENSG00000164104-HMGB2: "High mobility group box 2" and "Multifunctional protein with various roles in different cellular compartments. May act in a redox sensitive manner. In the nucleus is an abundant chromatin-associated non-histone protein involved in transcription, chromatin remodeling and V(D)J recombination and probably other processes. Binds DNA with a preference to non-canonical DNA structures such as single-stranded DNA. Can bent DNA and enhance DNA flexibility by looping thus providing a mechanism to promote activities on various gene promoters by enhancing transcription factor binding and/or bringing distant regulatory sequences into close proximity"
        'UBE2C', # correlated with TOP2A, MKI67, CENPA. https://www.proteinatlas.org/ENSG00000175063-UBE2C: "Ubiquitin conjugating enzyme E2 C" and "Acts as an essential factor of the anaphase promoting complex/cyclosome (APC/C), a cell cycle-regulated ubiquitin ligase that controls progression through mitosis. Acts by initiating 'Lys-11'-linked polyubiquitin chains on APC/C substrates, leading to the degradation of APC/C substrates by the proteasome and promoting mitotic exit."
        'KIF4A', # highly correlated with CENPE, CENPF. https://www.proteinatlas.org/ENSG00000090889-KIF4A: "Kinesin family member 4A" and "Iron-sulfur (Fe-S) cluster binding motor protein that has a role in chromosome segregation during mitosis 1. Translocates PRC1 to the plus ends of interdigitating spindle microtubules during the metaphase to anaphase transition, an essential step for the formation of an organized central spindle midzone and midbody and for successful cytokinesis 2, 3. May play a role in mitotic chromosomal positioning and bipolar spindle stabilization (By similarity)."
        'CCNA2', # highly correlated with TOP2A, MKI67. https://www.proteinatlas.org/ENSG00000145386-CCNA2: "Cyclin A2" and "Cyclin which controls both the G1/S and the G2/M transition phases of the cell cycle. Functions through the formation of specific serine/threonine protein kinase holoenzyme complexes with the cyclin-dependent protein kinases CDK1 or CDK2. The cyclin subunit confers the substrate specificity of these complexes and differentially interacts with and activates CDK1 and CDK2 throughout the cell cycle."
        'HMGA2', # https://www.proteinatlas.org/ENSG00000149948-HMGA2: "High mobility group AT-hook 2" and "Functions as a transcriptional regulator. Functions in cell cycle regulation through CCNA2. Plays an important role in chromosome condensation during the meiotic G2/M transition of spermatocytes". when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval, and high in 24h/36h/48h samples.
        'CCND2', # https://www.proteinatlas.org/ENSG00000118971-CCND2: "Cyclin D2" and "Regulatory component of the cyclin D2-CDK4 (DC) complex that phosphorylates and inhibits members of the retinoblastoma (RB) protein family including RB1 and regulates the cell-cycle during G(1)/S transition".
        'CCNB1', # https://www.proteinatlas.org/ENSG00000134057-CCNB1: "Cyclin B1" and "Essential for the control of the cell cycle at the G2/M (mitosis) transition"
        'CCNB2', # correlated with TOP2A, MKI67. https://www.proteinatlas.org/ENSG00000157456-CCNB2: "Cyclin B2" and "Essential for the control of the cell cycle at the G2/M (mitosis) transition."
        'NDC80', # https://www.proteinatlas.org/ENSG00000080986-NDC80: "NDC80 kinetochore complex component" and "Acts as a component of the essential kinetochore-associated NDC80 complex, which is required for chromosome segregation and spindle checkpoint activity"
        'DLGAP5', 
        'GTSE1', # highly correlated with TOP2A, MKI67. https://www.proteinatlas.org/ENSG00000075218-GTSE1: "G2 and S-phase expressed 1" and "The protein encoded by this gene is only expressed in the S and G2 phases of the cell cycle, where it colocalizes with cytoplasmic tubulin and microtubules. In response to DNA damage, the encoded protein accumulates in the nucleus and binds the tumor suppressor protein p53, shuttling it out of the nucleus and repressing its ability to induce apoptosis."
        'BIRC5', # correlated with MKI67. https://www.proteinatlas.org/ENSG00000089685-BIRC5: "Baculoviral IAP repeat containing 5" and "Multitasking protein that has dual roles in promoting cell proliferation and preventing apoptosis 1, 2, 3, 4, 5. Component of a chromosome passage protein complex (CPC) which is essential for chromosome alignment and segregation during mitosis and cytokinesis"
        'CDKN3', # correlated with MKI67. https://www.proteinatlas.org/ENSG00000100526-CDKN3: "Cyclin dependent kinase inhibitor 3" and "May play a role in cell cycle regulation"
        'PTTG1', # https://www.proteinatlas.org/ENSG00000164611-PTTG1: "PTTG1 regulator of sister chromatid separation, securin" and "The encoded protein is a homolog of yeast securin proteins, which prevent separins from promoting sister chromatid separation. It is an anaphase-promoting complex (APC) substrate that associates with a separin until activation of the APC. The gene product has transforming activity in vitro and tumorigenic activity in vivo, and the gene is highly expressed in various tumors" and "Regulatory protein, which plays a central role in chromosome stability, in the p53/TP53 pathway, and DNA repair".
        'PTTG1IP', # https://www.proteinatlas.org/ENSG00000183255-PTTG1IP: "PTTG1 interacting protein" and "May facilitate PTTG1 nuclear translocation".
        'POLE2', # https://www.proteinatlas.org/ENSG00000100479-POLE2: "DNA polymerase epsilon 2, accessory subunit" and "Accessory component of the DNA polymerase epsilon complex 1. Participates in DNA repair and in chromosomal DNA replication (By similarity)."
        'CLSPN', # correlated with PCNA. https://www.proteinatlas.org/ENSG00000092853-CLSPN: "Claspin" and "The product of this gene is an essential upstream regulator of checkpoint kinase 1 and triggers a checkpoint arrest of the cell cycle in response to replicative stress or DNA damage. The protein is also required for efficient DNA replication during a normal S phase."
        'PCLAF', # correlates with PCNA. https://www.proteinatlas.org/ENSG00000166803-PCLAF: "PCNA clamp associated factor" and "PCNA-binding protein that acts as a regulator of DNA repair during DNA replication. Following DNA damage, the interaction with PCNA is disrupted, facilitating the interaction between monoubiquitinated PCNA and the translesion DNA synthesis DNA polymerase eta (POLH) at stalled replisomes, facilitating the bypass of replication-fork-blocking lesions. Also acts as a regulator of centrosome number."
        'UHRF1', # correlates with MCM4. https://www.proteinatlas.org/ENSG00000276043-UHRF1: "Ubiquitin like with PHD and ring finger domains 1" and "Multidomain protein that acts as a key epigenetic regulator by bridging DNA methylation and chromatin modification. Specifically recognizes and binds hemimethylated DNA at replication forks via its YDG domain and recruits DNMT1 methyltransferase to ensure faithful propagation of the DNA methylation patterns through DNA replication."
        'CDC6', # correlates with H2AFZ, MCM10, MCM4. https://www.proteinatlas.org/ENSG00000094804-CDC6: "Cell division cycle 6" and "Involved in the initiation of DNA replication. Also participates in checkpoint controls that ensure DNA replication is completed before mitosis is initiated."
        'CCNE2', # https://www.proteinatlas.org/ENSG00000175305-CCNE2: "Cyclin E2" and "Essential for the control of the cell cycle at the late G1 and early S phase"
        'ZWINT', # https://www.proteinatlas.org/ENSG00000122952-ZWINT: "ZW10 interacting kinetochore protein" and "Part of the MIS12 complex, which is required for kinetochore formation and spindle checkpoint activity. Required to target ZW10 to the kinetochore at prometaphase"
        'DTL', # correlated with PCLAF, CENPU, MCM4, CLSPN. https://www.proteinatlas.org/ENSG00000143476-DTL: "Denticleless E3 ubiquitin protein ligase homolog" and "Substrate-specific adapter of a DCX (DDB1-CUL4-X-box) E3 ubiquitin-protein ligase complex required for cell cycle control, DNA damage response and translesion DNA synthesis. The DCX(DTL) complex, also named CRL4(CDT2) complex, mediates the polyubiquitination and subsequent degradation of CDT1, CDKN1A/p21(CIP1), FBH1, KMT5A and SDE2 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13. CDT1 degradation in response to DNA damage is necessary to ensure proper cell cycle regulation of DNA replication 14, 15, 16. CDKN1A/p21(CIP1) degradation during S phase or following UV irradiation is essential to control replication licensing"
        'HMMR', 
        'CEP70', # https://www.proteinatlas.org/ENSG00000114107-CEP70: "Centrosomal protein 70" and "Plays a role in the organization of both preexisting and nascent microtubules in interphase cells. During mitosis, required for the organization and orientation of the mitotic spindle."
        'CEP55', # correlated with MKI67. https://www.proteinatlas.org/ENSG00000138180-CEP55: "Centrosomal protein 55" and "Plays a role in mitotic exit and cytokinesis 1, 2. Recruits PDCD6IP and TSG101 to midbody during cytokinesis. Required for successful completion of cytokinesis"
        'CDK1', # correlated with TOP2A, NUSAP1, MKI67. https://www.proteinatlas.org/ENSG00000170312-CDK1: "Cyclin dependent kinase 1" and "Plays a key role in the control of the eukaryotic cell cycle by modulating the centrosome cycle as well as mitotic onset; promotes G2-M transition, and regulates G1 progress and G1-S transition via association with multiple interphase cyclins".
        'NCAPG', # correlated with MKI67 and TOP2A. https://www.proteinatlas.org/ENSG00000109805-NCAPG: "Non-SMC condensin I complex subunit G" and "Regulatory subunit of the condensin complex, a complex required for conversion of interphase chromatin into mitotic-like condense chromosomes. The condensin complex probably introduces positive supercoils into relaxed DNA in the presence of type I topoisomerases and converts nicked DNA into positive knotted forms in the presence of type II topoisomerases."
        'KIF14', # https://www.proteinatlas.org/ENSG00000118193-KIF14: "Kinesin family member 14" and "Microtubule motor protein that binds to microtubules with high affinity through each tubulin heterodimer and has an ATPase activity (By similarity). Plays a role in many processes like cell division, cytokinesis and also in cell proliferation and apoptosis 1, 2. During cytokinesis, targets to central spindle and midbody through its interaction with PRC1 and CIT respectively 3. Regulates cell growth through regulation of cell cycle progression and cytokinesis"
        'NASP', # https://www.proteinatlas.org/ENSG00000132780-NASP: "Required for DNA replication, normal cell cycle progression and cell proliferation. Forms a cytoplasmic complex with HSP90 and H1 linker histones and stimulates HSP90 ATPase activity. NASP and H1 histone are subsequently released from the complex and translocate to the nucleus where the histone is released for binding to DNA." and "This gene encodes a H1 histone binding protein that is involved in transporting histones into the nucleus of dividing cells. Multiple isoforms are encoded by transcript variants of this gene. The somatic form is expressed in all mitotic cells, is localized to the nucleus, and is coupled to the cell cycle."
        'TPX2', # https://www.proteinatlas.org/ENSG00000088325-TPX2: "TPX2 microtubule nucleation factor" and "Spindle assembly factor required for normal assembly of mitotic spindles. Required for normal assembly of microtubules during apoptosis. Required for chromatin and/or kinetochore dependent microtubule nucleation."

        'TLK1', # https://www.proteinatlas.org/ENSG00000198586-TLK1: "Tousled like kinase 1" and "Rapidly and transiently inhibited by phosphorylation following the generation of DNA double-stranded breaks during S-phase. This is cell cycle checkpoint and ATM-pathway dependent and appears to regulate processes involved in chromatin assembly".

        # nucleotide etc biosynthesis related, i think.
        'TYMS', # highly correlated with CLSPN and H2AFZ. https://www.proteinatlas.org/ENSG00000176890-TYMS: "Thymidylate synthetase" and "Thymidylate synthase catalyzes the methylation of deoxyuridylate to deoxythymidylate using, 10-methylenetetrahydrofolate (methylene-THF) as a cofactor. This function maintains the dTMP (thymidine-5-prime monophosphate) pool critical for DNA replication and repair."
        'TK1', # correlated with TYMS. https://www.proteinatlas.org/ENSG00000167900-TK1: "Thymidine kinase 1" and "The protein encoded by this gene is a cytosolic enzyme that catalyzes the addition of a gamma-phosphate group to thymidine. This creates dTMP and is the first step in the biosynthesis of dTTP, which is one component required for DNA replication. The encoded protein, whose levels fluctuate depending on the cell cycle stage, can act as a low activity dimer or a high activity tetramer."
        'GINS2', # correlated with TYMS. https://www.proteinatlas.org/ENSG00000131153-GINS2: "GINS complex subunit 2" and "The GINS complex plays an essential role in the initiation of DNA replication, and progression of DNA replication forks. GINS complex seems to bind preferentially to single-stranded DNA".
        'GINS4', # https://www.proteinatlas.org/ENSG00000147536-GINS4: "GINS complex subunit 4" and "The GINS complex plays an essential role in the initiation of DNA replication, and progression of DNA replication forks. GINS4 is important for GINS complex assembly. GINS complex seems to bind preferentially to single-stranded DNA".

        'CHEK1', # https://www.proteinatlas.org/ENSG00000149554-CHEK1: "Checkpoint kinase 1" and "Serine/threonine-protein kinase which is required for checkpoint-mediated cell cycle arrest and activation of DNA repair in response to the presence of DNA damage or unreplicated DNA 1, 2, 3, 4, 5, 6, 7, 8, 9. May also negatively regulate cell cycle progression during unperturbed cell cycles"
        'ATAD5', # https://www.proteinatlas.org/ENSG00000176208-ATAD5: "ATPase family AAA domain containing 5" and "Has an imporant role in DNA replication and in maintaining genome integrity during replication stress"

        # maybe more DNA repair related than cell cycle related? but still cell cycle related, i think.
        'RAD51AP1', # correlates with PCNA, CLSPN, H2AFZ, MCM10. https://www.proteinatlas.org/ENSG00000111247-RAD51AP1: "RAD51 associated protein 1" and "Structure-specific DNA-binding protein involved in DNA repair by promoting RAD51-mediated homologous recombination" and "Tissue cell type classificationi Cell type enriched (Colon - Mitotic cells (Colon), Lung - Mitotic cells (Lung), Skin - Mitotic cells (Skin), Stomach - Mitotic cells (Stomach), Testis - Early spermatids)"
        'XRCC2', # correlates with CLSPN in MCM4, H2AFZ. https://www.proteinatlas.org/ENSG00000196584-XRCC2: "X-ray repair cross complementing 2" and "Tissue cell type classificationi Cell type enriched (Lung - Mitotic cells (Lung), Skin - Mitotic cells (Skin), Stomach - Mitotic cells (Stomach), Testis - Early spermatids)"
        'BRCA1', # correlates with PCNA, TYMS, CLSPN, H2AFZ, MCM4. https://www.proteinatlas.org/ENSG00000012048-BRCA1: "BRCA1 DNA repair associated" and "E3 ubiquitin-protein ligase that specifically mediates the formation of 'Lys-6'-linked polyubiquitin chains and plays a central role in DNA repair by facilitating cellular responses to DNA damage" and "Regulates centrosomal microtubule nucleation 12. Required for appropriate cell cycle arrests after ionizing irradiation in both the S-phase and the G2 phase of the cell cycle"

        # 221210: SMC4 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'SMC4', # https://www.proteinatlas.org/ENSG00000113810-SMC4: "Structural maintenance of chromosomes 4" and "Central component of the condensin complex, a complex required for conversion of interphase chromatin into mitotic-like condense chromosomes. The condensin complex probably introduces positive supercoils into relaxed DNA in the presence of type I topoisomerases and converts nicked DNA into positive knotted forms in the presence of type II topoisomerases."

        # 221210: SMC4 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'SMC2', # https://www.proteinatlas.org/ENSG00000136824-SMC2: "Structural maintenance of chromosomes 2" and "Central component of the condensin complex, a complex required for conversion of interphase chromatin into mitotic-like condense chromosomes. The condensin complex probably introduces positive supercoils into relaxed DNA in the presence of type I topoisomerases and converts nicked DNA into positive knotted forms in the presence of type II topoisomerases."

        # 221210: each of MCM3, MCM6, MCM4, MCM2, MCM5 looked quite/very salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'MCM3', # https://www.proteinatlas.org/ENSG00000112118-MCM3: "Minichromosome maintenance complex component 3" and "Acts as component of the MCM2-7 complex (MCM complex) which is the putative replicative helicase essential for 'once per cell cycle' DNA replication initiation and elongation in eukaryotic cells."
        'MCM6', 'MCM4', 'MCM7', 'MCM2',
        'MCM5', # weirdly, MCM5 appeared far from the other MCMs above in the similarity matrix ordered according to the dendogram...

        'MCMBP', # https://www.proteinatlas.org/ENSG00000197771-MCMBP: "Minichromosome maintenance complex binding protein" and "Associated component of the MCM complex that acts as a regulator of DNA replication. Binds to the MCM complex during late S phase and promotes the disassembly of the MCM complex from chromatin, thereby acting as a key regulator of pre-replication complex (pre-RC) unloading from replicated DNA".

        'RGCC', # https://www.proteinatlas.org/ENSG00000102760-RGCC: "Regulator of cell cycle" and "This gene is thought to regulate cell cycle progression. It is induced by p53 in response to DNA damage, or by sublytic levels of complement system proteins that result in activation of the cell cycle. The encoded protein localizes to the cytoplasm during interphase and to centrosomes during mitosis. The protein forms a complex with polo-like kinase 1. The protein also translocates to the nucleus in response to treatment with complement system proteins, and can associate with and increase the kinase activity of cell division cycle 2 protein. In different assays and cell types, overexpression of this protein has been shown to activate or suppress cell cycle progression."

        
        'NUSAP1', # https://www.proteinatlas.org/ENSG00000137804-NUSAP1: "Nucleolar and spindle associated protein 1" and "Microtubule-associated protein with the capacity to bundle and stabilize microtubules (By similarity). May associate with chromosomes and promote the organization of mitotic spindle microtubules around them."

        # 'TOP2B', # https://www.proteinatlas.org/ENSG00000077097-TOP2B: "DNA topoisomerase II beta" and "This nuclear enzyme is involved in processes such as chromosome condensation, chromatid separation, and the relief of torsional stress that occurs during DNA transcription and replication." and "Plays a role in B-cell differentiation"
        # # 230227: the following were in a relatively strong module containing TOP2B (ordered by sub-modules - first of each was SPN, TOP2B, DKM1A and CHST12):
        # 'SPN', 'MAP7', 'CDK6',
        # # 'TOP2B',
        # 'CAPRIN1', 'E2F3',
        # 'FADS1', # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6949104/: "Suppression of FADS1 induces ROS generation, cell cycle arrest, and apoptosis in melanocytes"
        # 'KDM1A', 'EXTL2', 'CEP57L1', 'ZNF711', 'CCND2',
        # 'CHST12', 'FAM117A', 'TEC', 'MED12L', 'IGSF10', 'TRIM58', 'CEP70',

        'ANP32B', # https://www.proteinatlas.org/ENSG00000136938-ANP32B: "Acidic nuclear phosphoprotein 32 family member B" and "Multifunctional protein that is involved in the regulation of many processes including cell proliferation, apoptosis, cell cycle progression or transcription 1, 2. Regulates the proliferation of neuronal stem cells, differentiation of leukemic cells and progression from G1 to S phase of the cell cycle"

        'BTG1', # https://www.proteinatlas.org/ENSG00000133639-BTG1: "BTG anti-proliferation factor 1" and "This gene is a member of an anti-proliferative gene family that regulates cell growth and differentiation. Expression of this gene is highest in the G0/G1 phases of the cell cycle and downregulated when cells progressed through G1"
        'BTG2', # https://www.proteinatlas.org/ENSG00000159388-BTG2: "BTG anti-proliferation factor 2" and "This encoded protein is involved in the regulation of the G1/S transition of the cell cycle"
        'BTG3', # https://www.proteinatlas.org/ENSG00000154640-BTG3: "BTG anti-proliferation factor 3" and "Overexpression impairs serum-induced cell cycle progression from the G0/G1 to S phase". also, was upregulated in cells that waited for 24/36/48h.



        # 230227: the following were in a module containing PCNA (ordered by sub-modules):
        *['EBNA1BP2', 'GINS2', 'SLC29A1', 'CTNNAL1', 'MND1', 'NUSAP1'],
        *['UHRF1', 'MCM3', 'HELLS'],
        *['H2AFZ', 'RANBP1'],
        *['BRCA1', 'CLSPN', 'CDC6', 'MCM10', 'ZWINT', 'FEN1', 'PCNA', 'CENPK', 'SPC25', 'CDC45'],
        *['YBX1', 'CCT5', 'DTL', 'POLE2'],

        'CDK6', # https://www.proteinatlas.org/ENSG00000105810-CDK6: "Cyclin dependent kinase 6" and "Serine/threonine-protein kinase involved in the control of the cell cycle and differentiation; promotes G1/S transition. Phosphorylates pRB/RB1 and NPM1. Interacts with D-type G1 cyclins during interphase at G1 to form a pRB/RB1 kinase and controls the entrance into the cell cycle. Involved in initiation and maintenance of cell cycle exit during cell differentiation; prevents cell proliferation and regulates negatively cell differentiation, but is required for the proliferation of specific cell types (e.g. erythroid and hematopoietic cells)." NOTE: CDK6 Levels Regulate Quiescence Exit in Human Hematopoietic Stem Cells (https://www.cell.com/cell-stem-cell/pdf/S1934-5909(15)00018-1.pdf, 2015)

        # 'CDK13', #?

        'KMT2E', 

        'BOD1L1', # https://www.proteinatlas.org/ENSG00000038219-BOD1L1: "Biorientation of chromosomes in cell division 1 like 1" and "Component of the fork protection machinery required to protect stalled/damaged replication forks from uncontrolled DNA2-dependent resection. Acts by stabilizing RAD51 at stalled replication forks and protecting RAD51 nucleofilaments from the antirecombinogenic activities of FBH1 and BLM"
        'RAD21', # https://www.proteinatlas.org/ENSG00000164754-RAD21: "RAD21 cohesin complex component" and "[Double-strand-break repair protein rad21 homolog]: As a member of the cohesin complex, involved in sister chromatid cohesion from the time of DNA replication in S phase to their segregation in mitosis, a function that is essential for proper chromosome segregation, post-replicative DNA repair, and the prevention of inappropriate recombination between repetitive regions"
        'BAZ1A', # https://www.proteinatlas.org/ENSG00000198604-BAZ1A: "Bromodomain adjacent to zinc finger domain 1A" and "Regulatory subunit of the ATP-dependent ACF-1 and ACF-5 ISWI chromatin remodeling complexes, which form ordered nucleosome arrays on chromatin and slide edge- and center-positioned histone octamers away from their original location on the DNA template to facilitate access to DNA during DNA-templated processes such as DNA replication, transcription, and repair"
        'PAK2', # https://www.proteinatlas.org/ENSG00000180370-PAK2: "P21 (RAC1) activated kinase 2" and "Serine/threonine protein kinase that plays a role in a variety of different signaling pathways including cytoskeleton regulation, cell motility, cell cycle progression, apoptosis or proliferation.".

        'TPR', # https://www.proteinatlas.org/ENSG00000047410-TPR: "Translocated promoter region, nuclear basket protein" and "Biological process (UniProt)i Cell cycle, Cell division, Mitosis, mRNA transport, Protein transport, Translocation, Transport" and "Implicated in nuclear export of mRNAs transcribed from heat shock gene promoters; associates both with chromatin in the HSP70 promoter and with mRNAs transcribed from this promoter under stress-induced conditions" and "Functions as a scaffolding element in the nuclear phase of the NPC essential for normal nucleocytoplasmic transport of proteins and mRNAs, plays a role in the establishment of nuclear-peripheral chromatin compartmentalization in interphase, and in the mitotic spindle checkpoint signaling during mitosis."

        'MMS22L', # https://www.proteinatlas.org/ENSG00000146263-MMS22L: "MMS22 like, DNA repair protein" and "Component of the MMS22L-TONSL complex, a complex that promotes homologous recombination-mediated repair of double-strand breaks (DSBs) at stalled or collapsed replication forks 1, 2, 3, 4, 5, 6, 7. The MMS22L-TONSL complex is required to maintain genome integrity during DNA replication".

        'CHAF1A', # https://www.proteinatlas.org/ENSG00000167670-CHAF1A: "Chromatin assembly factor 1 subunit A" and "Core component of the CAF-1 complex, a complex that is thought to mediate chromatin assembly in DNA replication and DNA repair. Assembles histone octamers onto replicating DNA in vitro. CAF-1 performs the first step of the nucleosome assembly process, bringing newly synthesized histones H3 and H4 to replicating DNA; histones H2A/H2B can bind to this chromatin precursor subsequent to DNA replication to complete the histone octamer. It may play a role in heterochromatin maintenance in proliferating cells by bringing newly synthesized cbx proteins to heterochromatic DNA replication foci"

        'TIMELESS', # https://www.proteinatlas.org/ENSG00000111602-TIMELESS: "Timeless circadian regulator" and "Plays an important role in the control of DNA replication, maintenance of replication fork stability, maintenance of genome stability throughout normal DNA replication, DNA repair and in the regulation of the circadian clock 1, 2, 3, 4, 5. Required to stabilize replication forks during DNA replication by forming a complex with TIPIN"

        'VRK1', # https://www.proteinatlas.org/ENSG00000100749-VRK1: "VRK serine/threonine kinase 1" and "Serine/threonine kinase involved in cell cycle, nuclear condensation and transcription regulation 1, 2, 3. Involved in Golgi disassembly during the cell cycle: following phosphorylation by PLK3 during mitosis, required to induce Golgi fragmentation".

        'BRCA2', # https://www.proteinatlas.org/ENSG00000139618-BRCA2: "BRCA2 DNA repair associated" and "Involved in double-strand break repair and/or homologous recombination. Binds RAD51 and potentiates recombinational DNA repair by promoting assembly of RAD51 onto single-stranded DNA (ssDNA). Acts by targeting RAD51 to ssDNA over double-stranded DNA, enabling RAD51 to displace replication protein-A (RPA) from ssDNA and stabilizing RAD51-ssDNA filaments by blocking ATP hydrolysis" and "Binds selectively to ssDNA, and to ssDNA in tailed duplexes and replication fork structures. May play a role in the extension step after strand invasion at replication-dependent DNA double-strand breaks; together with PALB2 is involved in both POLH localization at collapsed replication forks and DNA polymerization activity. In concert with NPM1, regulates centrosome duplication"

        'MSH6', # https://www.proteinatlas.org/ENSG00000116062-MSH6: "MutS homolog 6" and "Component of the post-replicative DNA mismatch repair system (MMR). Heterodimerizes with MSH2 to form MutS alpha, which binds to DNA mismatches thereby initiating DNA repair. When bound, MutS alpha bends the DNA helix and shields approximately 20 base pairs, and recognizes single base mismatches and dinucleotide insertion-deletion loops (IDL) in the DNA"

        'RRM1', # https://www.proteinatlas.org/ENSG00000167325-RRM1: "Ribonucleotide reductase catalytic subunit M1" and "Provides the precursors necessary for DNA synthesis. Catalyzes the biosynthesis of deoxyribonucleotides from the corresponding ribonucleotides".

        'DUT', # https://www.proteinatlas.org/ENSG00000128951-DUT: "Deoxyuridine triphosphatase" and "hydrolyzes dUTP to dUMP and pyrophosphate. This reaction serves two cellular purposes: providing a precursor (dUMP) for the synthesis of thymine nucleotides needed for DNA replication, and limiting intracellular pools of dUTP. Elevated levels of dUTP lead to increased incorporation of uracil into DNA, which induces extensive excision repair mediated by uracil glycosylase. This repair process, resulting in the removal and reincorporation of dUTP, is self-defeating and leads to DNA fragmentation and cell death".

        'RPA3', # https://www.proteinatlas.org/ENSG00000106399-RPA3: "Replication protein A3" and "As part of the heterotrimeric replication protein A complex (RPA/RP-A), binds and stabilizes single-stranded DNA intermediates that form during DNA replication or upon DNA stress. It prevents their reannealing and in parallel, recruits and activates different proteins and complexes involved in DNA metabolism"

        'DNMT1', # https://www.proteinatlas.org/ENSG00000130816-DNMT1: "DNA methyltransferase 1" and "Methylates CpG residues. Preferentially methylates hemimethylated DNA. Associates with DNA replication sites in S phase maintaining the methylation pattern in the newly synthesized strand, that is essential for epigenetic inheritance. Associates with chromatin during G2 and M phases to maintain DNA methylation independently of replication"

        'CDT1', # https://www.proteinatlas.org/ENSG00000167513-CDT1: "Chromatin licensing and DNA replication factor 1" and "Required for both DNA replication and mitosis 1, 2, 3, 4, 5. DNA replication licensing factor, required for pre-replication complex assembly".

        'PRIM1', # https://www.proteinatlas.org/ENSG00000198056-PRIM1: "DNA primase subunit 1" and "Catalytic subunit of the DNA primase complex and component of the DNA polymerase alpha complex (also known as the alpha DNA polymerase-primase complex - primosome/replisome) which play an essential role in the initiation of DNA synthesis"

        'GMNN', # https://www.proteinatlas.org/ENSG00000112312-GMNN: "Geminin DNA replication inhibitor" and "Inhibits DNA replication by preventing the incorporation of MCM complex into pre-replication complex (pre-RC) 1, 2, 3, 4. It is degraded during the mitotic phase of the cell cycle 5, 6, 7. Its destruction at the metaphase-anaphase transition permits replication in the succeeding cell cycle"

        'NCAPG2', # https://www.proteinatlas.org/ENSG00000146918-NCAPG2: "Non-SMC condensin II complex subunit G2" and "Regulatory subunit of the condensin-2 complex, a complex which establishes mitotic chromosome architecture and is involved in physical rigidity of the chromatid axis"

        'RAN', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, seemed somewhat batchy, and clustered nicely with HMGB2 (among others). also, its function sounds lateral: https://www.proteinatlas.org/ENSG00000132341-RAN: "RAN, member RAS oncogene family" and "RAN (GTP-bound form) triggers microtubule assembly at mitotic chromosomes and is required for normal mitotic spindle assembly and chromosome segregation 26, 27. Required for normal progress through mitosis 28, 29, 30. The complex with BIRC5/survivin plays a role in mitotic spindle formation by serving as a physical scaffold to help deliver the RAN effector molecule TPX2 to microtubules" and "Biological process (UniProt)i: Cell cycle, Cell division, Host-virus interaction, Mitosis, Protein transport, Transport"

        'MAD2L1', # https://www.proteinatlas.org/ENSG00000164109-MAD2L1: "Mitotic arrest deficient 2 like 1" and "Component of the spindle-assembly checkpoint that prevents the onset of anaphase until all chromosomes are properly aligned at the metaphase plate"

        'ATAD2', # https://www.proteinatlas.org/ENSG00000156802-ATAD2: "ATPase family AAA domain containing 2" and "Single cell type expression clusteri: Non-specific - Cell proliferation (mainly)" and "Tissue cell type classificationi: Cell type enriched (Skin - Mitotic cells (Skin), Testis - Spermatocytes)" and "Involved in the estrogen-induced cell proliferation and cell cycle progression of breast cancer cells". also very highly correlated with MKI67 in a subset of B cells.
        'ARHGAP11A', # https://www.proteinatlas.org/ENSG00000198826-ARHGAP11A: "Rho GTPase activating protein 11A" and "Cell type enriched (Colon - Mitotic cells (Colon), Lung - Mitotic cells (Lung), Skin - Mitotic cells (Skin), Stomach - Mitotic cells (Stomach), Testis - Spermatocytes)". also very highly correlated with MKI67 in a subset of B cells.

        'KIF20B', # https://www.proteinatlas.org/ENSG00000138182-KIF20B: "Kinesin family member 20B" and "Plus-end-directed motor enzyme that is required for completion of cytokinesis 1, 2. Required for proper midbody organization and abscission in polarized cortical stem cells".
        'KIF11', # https://www.proteinatlas.org/ENSG00000138160-KIF11: "Kinesin family member 11" and "Motor protein required for establishing a bipolar spindle during mitosis 1. Required in non-mitotic cells for transport of secretory proteins from the Golgi complex to the cell surface".
        'KIF15', # https://www.proteinatlas.org/ENSG00000163808-KIF15: "Kinesin family member 15" and "Plus-end directed kinesin-like motor enzyme involved in mitotic spindle assembly".
        'KNL1', # https://www.proteinatlas.org/ENSG00000137812-KNL1: "Kinetochore scaffold 1" and "Performs two crucial functions during mitosis: it is essential for spindle-assembly checkpoint signaling and for correct chromosome alignment. Required for attachment of the kinetochores to the spindle microtubules".
        'BUB1B', # https://www.proteinatlas.org/ENSG00000156970-BUB1B: "BUB1 mitotic checkpoint serine/threonine kinase B" and "Essential component of the mitotic checkpoint. Required for normal mitosis progression. The mitotic checkpoint delays anaphase until all chromosomes are properly attached to the mitotic spindle"
        
        *REPLICATION_DEPENDENT_HISTONE_GENE_NAMES,
    ],

    'ribosomal': RIBOSOMAL_PROTEIN_GENE_NAMES,
    # 'highly_correlated_with_ribosomal': [
    'ribosomal_related': [
        'EEF1B2', # was very highly correlated with multiple RPL and RPS genes. https://www.proteinatlas.org/ENSG00000114942-EEF1B2: "Eukaryotic translation elongation factor 1 beta 2"
        'PRSS57', # was highly correlated with RPS24, RPLP0, RPL7. https://www.proteinatlas.org/ENSG00000185198-PRSS57: "Serine protease 57"
        'RACK1', # very highly correlated with RPS3A (0.97) and other ribosomal proteins. https://www.proteinatlas.org/ENSG00000204628-RACK1: "Receptor for activated C kinase 1" and "Component of the 40S ribosomal subunit involved in translational repression"
        'SNHG8', # highly correlated with RPS24 (0.88) and other ribosomal proteins. https://www.genecards.org/cgi-bin/carddisp.pl?gene=SNHG8: "Small Nucleolar RNA Host Gene 8"
        'HINT1', # was correlated (by eye, but clear) with donor_fraction of each of the donors in demux_28_11_21_1, a batch with high ribosomal protein umis in general. also, highly correlated with ribosomal protein genes. https://www.proteinatlas.org/ENSG00000169567-HINT1: "Histidine triad nucleotide binding protein 1" and "This gene encodes a protein that hydrolyzes purine nucleotide phosphoramidates substrates, including AMP-morpholidate, AMP-N-alanine methyl ester, AMP-alpha-acetyl lysine methyl ester, and AMP-NH2. The encoded protein interacts with these substrates via a histidine triad motif. This gene is considered a tumor suppressor gene".
        'NPM1', # highly correlated with RPLP0 (0.92) and other ribosomal proteins. https://www.proteinatlas.org/ENSG00000181163-NPM1: "Nucleophosmin 1" and "Involved in diverse cellular processes such as ribosome biogenesis, centrosome duplication, protein chaperoning, histone assembly, cell proliferation, and regulation of tumor suppressors p53/TP53 and ARF. Binds ribosome presumably to drive ribosome nuclear export. Associated with nucleolar ribonucleoprotein structures and bind single-stranded nucleic acids. Acts as a chaperonin for the core histones H3, H2B and H4."

        'NCL', # https://www.proteinatlas.org/ENSG00000115053-NCL: "Nucleolin is the major nucleolar protein of growing eukaryotic cells. It is found associated with intranucleolar chromatin and pre-ribosomal particles. It induces chromatin decondensation by binding to histone H1. It is thought to play a role in pre-rRNA transcription and ribosome assembly. May play a role in the process of transcriptional elongation. Binds RNA oligonucleotides with 5'-UUAGGG-3' repeats more tightly than the telomeric single-stranded DNA 5'- TTAGGG-3' repeats"

        'HNRNPA1', # highly correlated with RPLP0 (0.88) and other ribosomal proteins. https://www.proteinatlas.org/ENSG00000135486-HNRNPA1: "Heterogeneous nuclear ribonucleoprotein A1" and "Involved in the packaging of pre-mRNA into hnRNP particles, transport of poly(A) mRNA from the nucleus to the cytoplasm and modulation of splice site selection 1. Plays a role in the splicing of pyruvate kinase PKM by binding repressively to sequences flanking PKM exon 9, inhibiting exon 9 inclusion and resulting in exon 10 inclusion and production of the PKM M2 isoform 2. Binds to the IRES and thereby inhibits the translation of the apoptosis protease activating factor APAF1 3. May bind to specific miRNA hairpins" and "Biological process (UniProt)i	Host-virus interaction, mRNA processing, mRNA splicing, mRNA transport, Transport"

        'ZFAS1', # highly correlated with RPS24 (0.84) and other ribosomal proteins. https://www.genecards.org/cgi-bin/carddisp.pl?gene=ZFAS1: "This gene represents a snoRNA host gene that produces a non-coding RNA. Increased expression or amplification of this locus is associated with cancer progression and metastasis. This transcript regulates expression of genes involved in differentiation. It may act a molecular sponge for microRNAs."

        'FAU', # was correlated (by eye, but clear) with donor_fraction of each of the donors in demux_28_11_21_1, a batch with high ribosomal protein umis in general. also, highly correlated with ribosomal protein genes. https://www.proteinatlas.org/ENSG00000149806-FAU: "FAU ubiquitin like and ribosomal protein S30 fusion" and "[Ubiquitin-like protein FUBI]: May have pro-apoptotic activity" and "This gene is the cellular homolog of the fox sequence in the Finkel-Biskis-Reilly murine sarcoma virus (FBR-MuSV). It encodes a fusion protein consisting of the ubiquitin-like protein fubi at the N terminus and ribosomal protein S30 at the C terminus. It has been proposed that the fusion protein is post-translationally processed to generate free fubi and free ribosomal protein S30".
        'NACA', # was correlated (by eye, but clear) with donor_fraction of each of the donors in demux_28_11_21_1, a batch with high ribosomal protein umis in general. also, highly correlated with ribosomal protein genes. https://www.proteinatlas.org/ENSG00000196531-NACA: "Nascent polypeptide associated complex subunit alpha" and "This gene encodes a protein that associates with basic transcription factor 3 (BTF3) to form the nascent polypeptide-associated complex (NAC). This complex binds to nascent proteins that lack a signal peptide motif as they emerge from the ribosome, blocking interaction with the signal recognition particle (SRP) and preventing mistranslocation to the endoplasmic reticulum"

        'UBA52', # was correlated (by eye, but clear) with donor_fraction of each of the donors in demux_28_11_21_1, a batch with high ribosomal protein umis in general. also, highly correlated with ribosomal protein genes. https://www.proteinatlas.org/ENSG00000221983-UBA52: "Ubiquitin A-52 residue ribosomal protein fusion product 1" and "This gene encodes a fusion protein consisting of ubiquitin at the N terminus and ribosomal protein L40 at the C terminus, a C-terminal extension protein (CEP)."

        'ARL6IP1', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, clustered nicely with multiple genes including RPL36A, RPS17, RPL17.
        
        'LINC01410', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, clustered somewhat nicely with multiple genes including RPL36A, RPS17, RPL17.
        'AL031777.3', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, clustered somewhat with multiple genes including RPL36A, RPS17, RPL17.
    ],

    'presumably_stress_or_apoptosis_related': [
        *IFN_MODULE_IN_HIGH_IFN_METACELLS,

        'TP53',
        'MDM2', # https://www.proteinatlas.org/ENSG00000135679-MDM2: "MDM2 proto-oncogene" and "E3 ubiquitin-protein ligase that mediates ubiquitination of p53/TP53, leading to its degradation by the proteasome. Inhibits p53/TP53- and p73/TP73-mediated cell cycle arrest and apoptosis by binding its transcriptional activation domain"
        'TP53BP2', # https://www.proteinatlas.org/ENSG00000143514-TP53BP2: "Tumor protein p53 binding protein 2" and "Regulator that plays a central role in regulation of apoptosis and cell growth via its interactions with proteins such as TP53".
        'TP53INP1', # https://www.proteinatlas.org/ENSG00000164938-TP53INP1: "Tumor protein p53 inducible nuclear protein 1" and "Antiproliferative and proapoptotic protein involved in cell stress response which acts as a dual regulator of transcription and autophagy". also higher in samples that waited 24h/36h/48h.

        'RIF1', # https://www.proteinatlas.org/ENSG00000080345-RIF1: "Replication timing regulatory factor 1" and "Key regulator of TP53BP1 that plays a key role in the repair of double-strand DNA breaks (DSBs) in response to DNA damage: acts by promoting non-homologous end joining (NHEJ)-mediated repair of DSBs 1, 2. In response to DNA damage, interacts with ATM-phosphorylated TP53BP1"

        # https://en.wikipedia.org/wiki/Poly_(ADP-ribose)_polymerase: "a family of proteins involved in a number of cellular processes such as DNA repair, genomic stability, and programmed cell death"
        'PARP1', 'PARP3', 'PARP15', 'PARP14', 'PARP8', 'PARP12', 'PARP10', 'PARP11', 'PARP4', 'PARP2', 'PARP16', 'PARP6',
        'PARP9', # https://www.proteinatlas.org/ENSG00000138496-PARP9: "Poly(ADP-ribose) polymerase family member 9" and "ADP-ribosyltransferase which, in association with E3 ligase DTX3L, plays a role in DNA damage repair and in immune responses including interferon-mediated antiviral defenses"

        'ADAR', # https://www.proteinatlas.org/ENSG00000160710-ADAR: "Adenosine deaminase RNA specific" and "Catalyzes the hydrolytic deamination of adenosine to inosine in double-stranded RNA (dsRNA) referred to as A-to-I RNA editing" and "genetic stability in the case of RNA virus genomes by changing sequences during viral RNA replication" and "Biological process (UniProt)i Antiviral defense, Immunity, Innate immunity, mRNA processing, RNA-mediated gene silencing"

        'SAMD9', # correlated with ISG15 (0.85). https://www.proteinatlas.org/ENSG00000205413-SAMD9: "Sterile alpha motif domain containing 9" and "Cell line expression clusteri Non-specific - Antiviral immune response (mainly)" and "May play a role in the inflammatory response to tissue injury and the control of extra-osseous calcification, acting as a downstream target of TNF-alpha signaling"

        'SAMD9L', # correlated with ISG15 (0.86). https://www.proteinatlas.org/ENSG00000177409-SAMD9L: "Sterile alpha motif domain containing 9 like" and "Cell line expression clusteri Non-specific - Antiviral immune response (mainly)" and "This gene encodes a cytoplasmic protein that acts as a tumor suppressor but also plays a key role in cell proliferation and the innate immune response to viral infection."

        'RABGAP1L', # the 5 genes it was most correlated with (considering only metacells of interest, i.e., non-disease and stem) were SAMD9L, XAF1, EIF2AK2, SP100, PARP14. https://www.proteinatlas.org/ENSG00000152061-RABGAP1L: "RAB GTPase activating protein 1 like" and "GTP-hydrolysis activating protein (GAP) for small GTPase RAB22A, converting active RAB22A-GTP to the inactive form RAB22A-GDP 1. Plays a role in endocytosis and intracellular protein transport".

        'DDX58', # https://www.proteinatlas.org/ENSG00000107201-DDX58: "DExD/H-box helicase 58" and "Innate immune receptor that senses cytoplasmic viral nucleic acids and activates a downstream signaling cascade leading to the production of type I interferons and pro-inflammatory cytokines"

        'DDX60L', # correlated with ISG15 (0.8). https://www.proteinatlas.org/ENSG00000181381-DDX60L: "DExD/H-box 60 like" and "Cell line expression clusteri Non-specific - Antiviral immune response (mainly)"

        'DDX60', # https://www.proteinatlas.org/ENSG00000137628-DDX60: "DExD/H-box helicase 60" and "Positively regulates DDX58/RIG-I- and IFIH1/MDA5-dependent type I interferon and interferon inducible gene expression in response to viral infection."

        'HERC5', # https://www.proteinatlas.org/ENSG00000138646-HERC5: "HECT and RLD domain containing E3 ubiquitin protein ligase 5" and "Major E3 ligase for ISG15 conjugation. Acts as a positive regulator of innate antiviral response in cells induced by interferon.".

        'STAT2', # https://www.proteinatlas.org/ENSG00000170581-STAT2: "Signal transducer and activator of transcription 2" and "Signal transducer and activator of transcription that mediates signaling by type I interferons (IFN-alpha and IFN-beta)".

        'PLSCR1', # correlated with ISG15 (0.84). https://www.proteinatlas.org/ENSG00000188313-PLSCR1: "Phospholipid scramblase 1" and "Catalyzes calcium-induced ATP-independent rapid bidirectional and non-specific movement of phospholipids (lipid scrambling or lipid flip-flop) between the inner and outer leaflet of the plasma membrane resulting in collapse of the phospholipid asymmetry which leads to phosphatidylserine externalization on the cell surface" and "May play a role in the antiviral response of interferon (IFN) by amplifying and enhancing the IFN response through increased expression of select subset of potent antiviral genes"

        'BST2', # https://www.proteinatlas.org/ENSG00000130303-BST2: "Bone marrow stromal cell antigen 2" and "IFN-induced antiviral host restriction factor which efficiently blocks the release of diverse mammalian enveloped viruses by directly tethering nascent virions to the membranes of infected cells. Acts as a direct physical tether, holding virions to the cell membrane and linking virions to each other. The tethered virions can be internalized by endocytosis and subsequently degraded or they can remain on the cell surface. In either case, their spread as cell-free virions is restricted 1, 2, 3, 4, 5, 6, 7, 8, 9. Its target viruses belong to diverse families"

        'STAT1', # https://www.proteinatlas.org/ENSG00000115415-STAT1: "Signal transducer and activator of transcription 1" and "Signal transducer and transcription activator that mediates cellular responses to interferons (IFNs), cytokine KITLG/SCF and other cytokines and other growth factors. Following type I IFN (IFN-alpha and IFN-beta) binding to cell surface receptors, signaling via protein kinases leads to activation of Jak kinases (TYK2 and JAK1) and to tyrosine phosphorylation of STAT1 and STAT2. The phosphorylated STATs dimerize and associate with ISGF3G/IRF-9 to form a complex termed ISGF3 transcription factor, that enters the nucleus 1. ISGF3 binds to the IFN stimulated response element (ISRE) to activate the transcription of IFN-stimulated genes (ISG), which drive the cell in an antiviral state. In response to type II IFN (IFN-gamma), STAT1 is tyrosine- and serine-phosphorylated 2. It then forms a homodimer termed IFN-gamma-activated factor (GAF), migrates into the nucleus and binds to the IFN gamma activated sequence (GAS) to drive the expression of the target genes, inducing a cellular antiviral state."

        'OAS2', # https://www.proteinatlas.org/ENSG00000111335-OAS2: "2'-5'-oligoadenylate synthetase 2" and "Interferon-induced, dsRNA-activated antiviral enzyme which plays a critical role in cellular innate antiviral response 1, 2. Activated by detection of double stranded RNA (dsRNA)"
        'OAS1', # https://www.proteinatlas.org/ENSG00000089127-OAS1: "2'-5'-oligoadenylate synthetase 1" and "Interferon-induced, dsRNA-activated antiviral enzyme which plays a critical role in cellular innate antiviral response 1. In addition, it may also play a role in other cellular processes such as apoptosis, cell growth, differentiation and gene regulation"
        'OAS3', # https://www.proteinatlas.org/ENSG00000111331-OAS3: "2'-5'-oligoadenylate synthetase 3" and "Interferon-induced, dsRNA-activated antiviral enzyme which plays a critical role in cellular innate antiviral response. In addition, it may also play a role in other cellular processes such as apoptosis, cell growth, differentiation and gene regulation."
        'OASL', # highly correlated with ISG15, IFIT3, IFIT1. https://www.proteinatlas.org/ENSG00000135114-OASL: "2'-5'-oligoadenylate synthetase like" and "Does not have 2'-5'-OAS activity, but can bind double-stranded RNA. Displays antiviral activity against encephalomyocarditis virus (EMCV) and hepatitis C virus (HCV) via an alternative antiviral pathway independent of RNase L."

        'MT2A', # very highly correlated with IFITM3 (0.93). https://www.proteinatlas.org/ENSG00000125148-MT2A: "Metallothionein 2A" and "This gene is a member of the metallothionein family of genes. Proteins encoded by this gene family are low in molecular weight, are cysteine-rich, lack aromatic residues, and bind divalent heavy metal ions, altering the intracellular concentration of heavy metals in the cell. These proteins act as anti-oxidants, protect against hydroxyl free radicals, are important in homeostatic control of metal in the cell, and play a role in detoxification of heavy metals. The encoded protein interacts with the protein encoded by the homeobox containing 1 gene in some cell types, controlling intracellular zinc levels, affecting apoptotic and autophagy pathways."

        'GBP1', # highly correlated with IFITM3 (0.86). https://www.proteinatlas.org/ENSG00000117228-GBP1: "Guanylate binding protein 1" and "Hydrolyzes GTP to GMP in 2 consecutive cleavage reactions (By similarity). Confers protection to several pathogens, including the bacterial pathogens Listeria monocytogenes and Mycobacterium bovis BCG as well as the protozoan pathogen Toxoplasma gondii (By similarity). Promotes IFN-gamma-mediated host defense against bacterial infections by regulating bacteriolytic peptide generation via its interaction with ubiquitin-binding protein SQSTM1, which delivers monoubiquitylated proteins to autolysosomes for the generation of bacteriolytic peptides (By similarity). Exhibits antiviral activity against influenza virus"

        'LAP3', # highly correlated with IFITM3 (0.87). https://www.proteinatlas.org/ENSG00000002549-LAP3: "Leucine aminopeptidase 3" and "Cell line expression clusteri Non-specific - Antiviral immune response (mainly)" and "Cytosolic metallopeptidase that catalyzes the removal of unsubstituted N-terminal hydrophobic amino acids from various peptides. The presence of Zn(2+) ions is essential for the peptidase activity, and the association with other cofactors can modulate the substrate spectificity of the enzyme. For instance, in the presence of Mn(2+), it displays a specific Cys-Gly hydrolyzing activity of Cys-Gly-S-conjugates. Involved in the metabolism of glutathione and in the degradation of glutathione S-conjugates, which may play a role in the control of the cell redox status."

        'LY6E', # highly correlated with ISG15 (0.88). https://www.proteinatlas.org/ENSG00000160932-LY6E: "Lymphocyte antigen 6 family member E" and "Cell line expression clusteri Non-specific - Antiviral immune response (mainly)" and "GPI-anchored cell surface protein that regulates T-lymphocytes proliferation, differentiation, and activation. Regulates the T-cell receptor (TCR) signaling by interacting with component CD3Z/CD247 at the plasma membrane, leading to CD3Z/CD247 phosphorylation modulation (By similarity). Restricts the entry of human coronaviruses, including SARS-CoV, MERS-CoV and SARS-CoV-2, by interfering with spike protein-mediated membrane fusion"

        'RNF213', # https://www.proteinatlas.org/ENSG00000173821-RNF213: "Ring finger protein 213" and "Atypical E3 ubiquitin ligase that can catalyze ubiquitination of both proteins and lipids, and which is involved in various processes, such as lipid metabolism, angiogenesis and cell-autonomous immunity 1, 2, 3, 4, 5, 6, 7. Acts as a key immune sensor by catalyzing ubiquitination of the lipid A moiety of bacterial lipopolysaccharide (LPS) via its RZ-type zinc-finger: restricts the proliferation of cytosolic bacteria, such as Salmonella, by generating the bacterial ubiquitin coat through the ubiquitination of LPS". https://www.nature.com/articles/s41467-021-26061-w: "RNF213 is a poorly characterized, interferon-induced megaprotein that is frequently mutated in Moyamoya disease, a rare cerebrovascular disorder. We report that interferon induces ISGylation and oligomerization of RNF213 on lipid droplets, where it acts as a sensor for ISGylated proteins. We show that RNF213 has broad antimicrobial activity in vitro and in vivo, counteracting infection with Listeria monocytogenes, herpes simplex virus 1, human respiratory syncytial virus and coxsackievirus B3, and we observe a striking co-localization of RNF213 with intracellular bacteria."

        'SP100', # https://www.proteinatlas.org/ENSG00000067066-SP100: "SP100 nuclear antigen" and "Together with PML, this tumor suppressor is a major constituent of the PML bodies, a subnuclear organelle involved in a large number of physiological processes including cell growth, differentiation and apoptosis.". https://pubmed.ncbi.nlm.nih.gov/8810287/: "Expression of the nuclear domain-associated proteins Sp100, PML, and NDP52, is enhanced by interferons (IFNs) on the mRNA and protein level"
        'SP110', # https://www.proteinatlas.org/ENSG00000135899-SP110: "SP110 nuclear body protein" and "Gene namei SP110 (IFI41, IFI75)" and "Cell line expression clusteri Non-specific - Antiviral immune response (mainly)" and "Transcription factor. May be a nuclear hormone receptor coactivator. Enhances transcription of genes with retinoic acid response elements (RARE)". https://jbiomedsci.biomedcentral.com/articles/10.1186/s12929-018-0434-4: "SP110, an interferon-induced nuclear protein, belongs to the SP100/SP140 protein family. Very recently, we showed that SP110b, an SP110 isoform, controls host innate immunity to Mycobacterium tuberculosis infection by regulating nuclear factor-κB (NF-κB) activity."

        'DUSP1', # https://www.proteinatlas.org/ENSG00000120129-DUSP1: "Dual specificity phosphatase 1" and "Biological process (UniProt)i Cell cycle, Stress response" and "This protein appears to play an important role in the human cellular response to environmental stress as well as in the negative regulation of cellular proliferation"

        'SLFN5', # https://www.proteinatlas.org/ENSG00000166750-SLFN5: "Schlafen family member 5" and "May have a role in hematopoietic cell differentiation". https://www.sciencedirect.com/science/article/pii/S0021925820605622: "Our results demonstrate that SLFN5, and to a lesser degree SLFN11 and SLFN12, are induced by IFNα in normal melanocytes, although in IFN-sensitive melanoma cells (SKMEL2, SKMEL5, and SKMEL28) only SLFN5 is induced."

        # 'PDCD4', # https://www.proteinatlas.org/ENSG00000150593-PDCD4: "Programmed cell death 4" and "Inhibits translation initiation and cap-dependent translation. May excert its function by hindering the interaction between EIF4A1 and EIF4G. Inhibits the helicase activity of EIF4A. Modulates the activation of JUN kinase. Down-regulates the expression of MAP4K1, thus inhibiting events important in driving invasion, namely, MAPK85 activation and consequent JUN-dependent transcription. May play a role in apoptosis".

        # 221210: IFI6 looked very salt and pepper in the enrichment view in mcview (interestingly, this is not the case for nimrodra's data, even though IFI6 was not marked as lateral. though maybe enough genes that are part of the interferon respons were marked as lateral, so that metacells didn't form more due to interferon response??). given the salt and pepper look and the relation to apoptosis, decided to mark as lateral.
        'IFI6', # https://www.proteinatlas.org/ENSG00000126709-IFI6: "Interferon alpha inducible protein 6" and "Plays a role in apoptosis, negatively regulating the intrinsinc apoptotic signaling pathway and TNFSF10-induced apoptosis 1, 2, 3. However, it has also been shown to have a pro-apoptotic activity 4. Has an antiviral activity towards hepatitis C virus/HCV by inhibiting the EGFR signaling pathway, which activation is required for entry of the virus into cells" and "may play a critical role in the regulation of apoptosis"

        # # 221210: TSPO correlated with IFI6, but https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0185767: "The macrophage marker translocator protein (TSPO)"
        # # 'TSPO',

        # 221210: IFIT1 looked very salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'IFIT1', # https://www.proteinatlas.org/ENSG00000185745-IFIT1: "Interferon induced protein with tetratricopeptide repeats 1" and "Interferon-induced antiviral RNA-binding protein that specifically binds single-stranded RNA bearing a 5'-triphosphate group (PPP-RNA), thereby acting as a sensor of viral single-stranded RNAs and inhibiting expression of viral messenger RNAs." and "originally identified as induced upon treatment with interferon"

        # 221210: IFIT3 looked very salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'IFIT3', # https://www.proteinatlas.org/ENSG00000119917-IFIT3: "Interferon induced protein with tetratricopeptide repeats 3" and "IFN-induced antiviral protein which acts as an inhibitor of cellular as well as viral processes, cell migration, proliferation, signaling, and viral replication."

        # EPSTI1 correlated with IFIT1 and IFIT3, but it seems specific to macrophages (https://www.proteinatlas.org/ENSG00000133106-EPSTI1). but monocytes are easy to identify, and anyway we are interested in HSPCs, so not a problem to mark it as lateral...
        'EPSTI1', # https://www.proteinatlas.org/ENSG00000133106-EPSTI1: "Epithelial stromal interaction 1" and "Plays a role in M1 macrophage polarization and is required for the proper regulation of gene expression during M1 versus M2 macrophage differentiation (By similarity). Might play a role in RELA/p65 and STAT1 phosphorylation and nuclear localization upon activation of macrophages (By similarity)."

        'EIF2AK2', # https://www.proteinatlas.org/ENSG00000055332-EIF2AK2: "Eukaryotic translation initiation factor 2 alpha kinase 2" and "IFN-induced dsRNA-dependent serine/threonine-protein kinase that phosphorylates the alpha subunit of eukaryotic translation initiation factor 2 (EIF2S1/eIF-2-alpha) and plays a key role in the innate immune response to viral infection"

        # decided to not mark IFI16 as lateral because https://www.proteinatlas.org/ENSG00000163565-IFI16 says "Could have a role in the regulation of hematopoietic differentiation". 230220: seems like there is no choice. it seems to be involved in interfon response - correlated with multiple interferon induced genes.
        'IFI16', # https://www.proteinatlas.org/ENSG00000163565-IFI16: "Could have a role in the regulation of hematopoietic differentiation through activation of unknown target genes. Controls cellular proliferation by modulating the functions of cell cycle regulatory factors including p53/TP53 and the retinoblastoma protein." and "This gene encodes a member of the HIN-200 (hematopoietic interferon-inducible nuclear antigens with 200 amino acid repeats) family of cytokines. The encoded protein contains domains involved in DNA binding, transcriptional regulation, and protein-protein interactions. The protein localizes to the nucleoplasm and nucleoli, and interacts with p53 and retinoblastoma-1. It modulates p53 function, and inhibits cell growth in the Ras/Raf signaling pathway."

        # 221210: IFI44L looked very salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'IFI44L', # https://www.proteinatlas.org/ENSG00000137959-IFI44L: "Interferon induced protein 44 like" and "Type I interferon-stimulated gene (ISG) that plays a critical role in antiviral and antibacterial activity"

        'IFI44', # https://www.proteinatlas.org/ENSG00000137965-IFI44: "Interferon induced protein 44" and "This protein aggregates to form microtubular structures"

        # 221210: MX1 looked very salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'MX1', # https://www.proteinatlas.org/ENSG00000157601-MX1: "Interferon-induced dynamin-like GTPase with antiviral activity against a wide range of RNA viruses and some DNA viruses"

        'TRIM22', # https://www.proteinatlas.org/ENSG00000132274-TRIM22: "Tripartite motif containing 22" and "Interferon-induced antiviral protein involved in cell innate immunity. The antiviral activity could in part be mediated by TRIM22-dependent ubiquitination of viral proteins. Plays a role in restricting the replication of HIV-1, encephalomyocarditis virus (EMCV) and hepatitis B virus (HBV)."

        # 221210: MX2 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'MX2', # https://www.proteinatlas.org/ENSG00000183486-MX2: "Interferon-induced dynamin-like GTPase with potent antiviral activity against human immunodeficiency virus type 1 (HIV-1). Acts by targeting the viral capsid and affects the nuclear uptake and/or stability of the HIV-1 replication complex and the subsequent chromosomal integration of the proviral DNA. Exhibits antiviral activity also against simian immunodeficiency virus (SIV-mnd). May play a role in regulating nucleocytoplasmic transport and cell-cycle progression."

        # 221210: IFIT2 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'IFIT2', # https://www.proteinatlas.org/ENSG00000119922-IFIT2: "Interferon induced protein with tetratricopeptide repeats 2" and "IFN-induced antiviral protein which inhibits expression of viral messenger RNAs lacking 2'-O-methylation of the 5' cap."

        # 221210: TRIM22 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'TRIM22', # https://www.proteinatlas.org/ENSG00000132274-TRIM22: "Tripartite motif containing 22" and "Interferon-induced antiviral protein involved in cell innate immunity"

        # 230116: TRIM73 looked very salt and pepper in mcview.
        'TRIM73', # https://www.proteinatlas.org/ENSG00000178809-TRIM73: "Tripartite motif containing 73". https://bmccancer.biomedcentral.com/articles/10.1186/s12885-022-09727-6 (2022): "Tripartite motif containing 73 (TRIM73) is a member of the TRIM family and is also called Tripartite motif-containing protein 50B (TRIM50B). It has been reported that the TRIM73 might act as an E3 ubiquitin ligase [8]. However, literatures about the TRIM73 are still rare, especially its biological function in tumors. Up to now, only Li et al. reported that the supermethylation of TRIM73 in plasma could be applied as an important indicator for the early diagnosis of pancreatic cancer [9]. Actually, the clinical values of TRIM family genes in tumorigenesis, development and prognosis have been manifested by a number of relevant studies and most of members of the TRIM family play a role in tumors as proto-oncogenes or tumor-promoting genes [10,11,12,13,14,15]."


        'TRIM56', # was correlated with ISG15 in Dissimilar-high-IFN metacells. https://www.proteinatlas.org/ENSG00000169871-TRIM56: "Tripartite motif containing 56" and "E3 ubiquitin-protein ligase that plays a key role in innate antiviral immunity by mediating ubiquitination of CGAS and STING1 1, 2. In response to pathogen- and host-derived double-stranded DNA (dsDNA), targets STING1 to 'Lys-63'-linked ubiquitination, thereby promoting its homodimerization, a step required for the production of type I interferon IFN-beta (By similarity)."

        'APOL6', # separated well between Dissimilar-high-IFN and normal metacells, and was correlated with ISG15 in Dissimilar-high-IFN metacells. https://www.proteinatlas.org/ENSG00000221963-APOL6: "Apolipoprotein L6". https://www.ncbi.nlm.nih.gov/gene/80830: "APOL6 has a high ranking for tumor aggressiveness, and interacts with genes in the STAT pathway, and with genes annotated as being part of the defense response, and the interferon signaling pathway."

        'TRIM38', # separated well between Dissimilar-high-IFN and normal metacells. https://www.proteinatlas.org/ENSG00000112343-TRIM38: "Tripartite motif containing 38" and "E3 ubiquitin-protein and E3 SUMO-protein ligase that acts as a regulator of innate immunity 1. Acts as a negative regulator of type I interferon IFN-beta production by catalyzing 'Lys-48'-linked polyubiquitination of AZI2/NAP1, leading to its degradation (By similarity). Mediates 'Lys-48'-linked polyubiquitination and proteasomal degradation of the critical TLR adapter TICAM1, inhibiting TLR3-mediated type I interferon signaling 2. Acts as positive regulator of the cGAS-STING pathway by acting as a E3 SUMO-protein ligase: mediates sumoylation of CGAS and STING, preventing their degradation and thereby activating the innate immune response to DNA virus (By similarity)".
        
        
        'MCL1', # separated well between Dissimilar-high-IFN and normal metacells. https://www.proteinatlas.org/ENSG00000143384-MCL1: "MCL1 apoptosis regulator, BCL2 family member" and "Involved in the regulation of apoptosis versus cell survival, and in the maintenance of viability but not of proliferation. Mediates its effects by interactions with a number of other regulators of apoptosis. Isoform 1 inhibits apoptosis. Isoform 2 promotes apoptosis".

        # 230116: MIF looked salt and pepper in mcview.
        'MIF', # https://www.proteinatlas.org/ENSG00000240972-MIF: "Macrophage migration inhibitory factor" and "Pro-inflammatory cytokine involved in the innate immune response to bacterial pathogens 1, 2, 3. The expression of MIF at sites of inflammation suggests a role as mediator in regulating the function of macrophages in host defense 4, 5, 6. Counteracts the anti-inflammatory activity of glucocorticoids"

        'IRF7', # https://www.proteinatlas.org/ENSG00000185507-IRF7: "Interferon regulatory factor 7" and "Key transcriptional regulator of type I interferon (IFN)- dependent immune responses and plays a critical role in the innate immune response against DNA and RNA viruses. Regulates the transcription of type I IFN genes (IFN-alpha and IFN-beta) and IFN-stimulated genes (ISG) by binding to an interferon-stimulated response element (ISRE) in their promoters 1, 2. Can efficiently activate both the IFN-beta (IFNB) and the IFN-alpha (IFNA) genes and mediate their induction via both the virus-activated, MyD88-independent pathway and the TLR-activated, MyD88-dependent pathway"

        'IRF1', # https://www.proteinatlas.org/ENSG00000125347-IRF1: "Interferon regulatory factor 1" and "Transcriptional regulator which displays a remarkable functional diversity in the regulation of cellular responses 1, 2, 3, 4, 5, 6, 7, 8, 9. Regulates transcription of IFN and IFN-inducible genes, host response to viral and bacterial infections, regulation of many genes expressed during hematopoiesis, inflammation, immune responses and cell proliferation and differentiation, regulation of the cell cycle and induction of growth arrest and programmed cell death following DNA damage" and "Its target genes for transcriptional activation activity include: genes involved in anti-viral response, such as IFN-alpha/beta, DDX58/RIG-I, TNFSF10/TRAIL, ZBP1, OAS1/2, PIAS1/GBP, EIF2AK2/PKR and RSAD2/viperin; antibacterial response, such as NOS2/INOS; anti-proliferative response, such as p53/TP53, LOX and CDKN1A; apoptosis, such as BBC3/PUMA, CASP1, CASP7 and CASP8; immune response, such as IL7, IL12A/B and IL15, PTGS2/COX2 and CYBB; DNA damage responses and DNA repair, such as POLQ/POLH; MHC class I expression, such as TAP1, PSMB9/LMP2, PSME1/PA28A, PSME2/PA28B and B2M and MHC class II expression, such as CIITA; metabolic enzymes, such as ACOD1/IRG1"

        'RIOK3', # https://www.proteinatlas.org/ENSG00000101782-RIOK3: "RIO kinase 3" and "Involved in regulation of type I interferon (IFN)-dependent immune response which plays a critical role in the innate immune response against DNA and RNA viruses".

        'PSMB9', # https://www.proteinatlas.org/ENSG00000240065-PSMB9: "Proteasome 20S subunit beta 9" and "This subunit is involved in antigen processing to generate class I binding peptides. Replacement of PSMB6 by PSMB9 increases the capacity of the immunoproteasome to cleave model peptides after hydrophobic and basic residues."
        'PSME2', # https://www.proteinatlas.org/ENSG00000100911-PSME2: "Proteasome activator subunit 2" and "Implicated in immunoproteasome assembly and required for efficient antigen processing. The PA28 activator complex enhances the generation of class I binding peptides by altering the cleavage pattern of the proteasome".

        # 221210: ISG15 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'ISG15', # https://www.genecards.org/cgi-bin/carddisp.pl?gene=ISG15: "Interferon-Stimulated Protein, 15 KDa". https://www.proteinatlas.org/ENSG00000187608-ISG15: "Ubiquitin-like protein which plays a key role in the innate immune response to viral infection either via its conjugation to a target protein (ISGylation) or via its action as a free or unconjugated protein. ISGylation involves a cascade of enzymatic reactions involving E1, E2, and E3 enzymes which catalyze the conjugation of ISG15 to a lysine residue in the target protein 1. Its target proteins include IFIT1, MX1/MxA, PPM1B, UBE2L6, UBA7, CHMP5, CHMP2A, CHMP4B and CHMP6. Isgylation of the viral sensor IFIH1/MDA5 promotes IFIH1/MDA5 oligomerization and triggers activation of innate immunity against a range of viruses, including coronaviruses, flaviviruses and picornaviruses"

        # 221210: ISG20 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'ISG20', # https://www.proteinatlas.org/ENSG00000172183-ISG20: "Interferon-induced antiviral exoribonuclease that acts on single-stranded RNA and also has minor activity towards single-stranded DNA. Exhibits antiviral activity against RNA viruses including hepatitis C virus (HCV), hepatitis A virus (HAV) and yellow fever virus (YFV) in an exonuclease-dependent manner. May also play additional roles in the maturation of snRNAs and rRNAs, and in ribosome biogenesis."

        # 221210: RSAD2 looked very salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'RSAD2', # https://www.proteinatlas.org/ENSG00000134321-RSAD2: "Radical S-adenosyl methionine domain containing 2" and "Interferon-inducible antiviral protein which plays a major role in the cell antiviral state induced by type I and type II interferon 1. Catalyzes the conversion of cytidine triphosphate (CTP) to 3'-deoxy-3',4'-didehydro-CTP (ddhCTP) via a SAM-dependent radical mechanism 2, 3. In turn, ddhCTP acts as a chain terminator for the RNA-dependent RNA polymerases from multiple viruses and directly inhibits viral replication 4. Therefore, inhibits a wide range of DNA and RNA viruses"

        'CLU', # https://www.proteinatlas.org/ENSG00000120885-CLU: "Clusterin" and "[Isoform 1]: Functions as extracellular chaperone that prevents aggregation of non native proteins 1, 2. Prevents stress-induced aggregation of blood plasma proteins" and "A mitochondrial form suppresses BAX-dependent release of cytochrome c into the cytoplasm and inhibit apoptosis 18, 19. Plays a role in the regulation of cell proliferation 20. An intracellular form suppresses stress-induced apoptosis by stabilizing mitochondrial membrane integrity through interaction with HSPA5" and "Plays a role in the clearance of immune complexes that arise during cell injury (By similarity)"

        'DTX3L', # https://www.proteinatlas.org/ENSG00000163840-DTX3L: "Deltex E3 ubiquitin ligase 3L" and "E3 ubiquitin-protein ligase which, in association with ADP-ribosyltransferase PARP9, plays a role in DNA damage repair and in interferon-mediated antiviral responses"

        'UBE2L6', # https://www.proteinatlas.org/ENSG00000156587-UBE2L6: "Ubiquitin conjugating enzyme E2 L6" and "Catalyzes the covalent attachment of ubiquitin or ISG15 to other proteins. Functions in the E6/E6-AP-induced ubiquitination of p53/TP53. Promotes ubiquitination and subsequent proteasomal degradation of FLT3".

        'TRIM25', # https://www.proteinatlas.org/ENSG00000121060-TRIM25: "Tripartite motif containing 25" and "Functions as a ubiquitin E3 ligase and as an ISG15 E3 ligase 1. Involved in innate immune defense against viruses by mediating ubiquitination of DDX58 and IFIH1"

        'CMPK2', # was highly correlated with IFIT1 (0.89), IFIT3 (0.89). https://www.proteinatlas.org/ENSG00000134326-CMPK2: "Cytidine/uridine monophosphate kinase 2" and "May participate in dUTP and dCTP synthesis in mitochondria".

        # 221210: JUN looked very salt and pepper in the enrichment view in mcview, and due to its presumable function (see below AP-1 transcription factor), decided to mark it as lateral.
        'JUN', # https://www.proteinatlas.org/ENSG00000177606-JUN: "Jun proto-oncogene, AP-1 transcription factor subunit" and https://www.nature.com/articles/ncb0502-e131: "The transcription factor AP-1 (activator protein-1) is involved in cellular proliferation, transformation and death. Using mice and cells lacking AP-1 components, the target-genes and molecular mechanisms mediating these processes were recently identified. Interestingly, the growth-promoting activity of c-Jun is mediated by repression of tumour suppressors, as well as upregulation of positive cell cycle regulators. Mostly, c-Jun is a positive regulator of cell proliferation, whereas JunB has the converse effect" (also related: https://www.nature.com/articles/nrc1209).

        # 221210: FOS looked very salt and pepper in the enrichment view in mcview, and due to its presumable function (see above AP-1 transcription factor), decided to mark it as lateral.
        'FOS', # https://www.proteinatlas.org/ENSG00000170345-FOS: "Fos proto-oncogene, AP-1 transcription factor subunit"
        'AHR', # correlated with FOS. https://www.proteinatlas.org/ENSG00000106546-AHR: "Aryl hydrocarbon receptor" and "Ligand-activated transcription factor that enables cells to adapt to changing conditions by sensing compounds from the environment, diet, microbiome and cellular metabolism"
        # 221210: FOSB looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function (see above AP-1 transcription factor), decided to mark it as lateral.
        'FOSB', # https://www.proteinatlas.org/ENSG00000125740-FOSB: "FosB proto-oncogene, AP-1 transcription factor subunit"

        # 221210: JUND looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function (see above AP-1 transcription factor), decided to mark it as lateral.
        'JUND', # https://www.proteinatlas.org/ENSG00000130522-JUND: "JunD proto-oncogene, AP-1 transcription factor subunit"

        # 221210: JUNB looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function (see above AP-1 transcription factor), decided to mark it as lateral.
        'JUNB', # https://www.proteinatlas.org/ENSG00000171223-JUNB: "JunB proto-oncogene, AP-1 transcription factor subunit"
        'IER2', # highly correlated with JUNB. https://www.proteinatlas.org/ENSG00000160888-IER2: "Immediate early response 2"

        'SORL1', # correlated with FOS. https://www.proteinatlas.org/ENSG00000137642-SORL1: "Sortilin related receptor 1" and "Sorting receptor that directs several proteins to their correct location within the cell (Probable). Along with AP-1 complex, involved Golgi apparatus - endosome sorting" and "Metabolic regulator, which functions to maintain the adequate balance between lipid storage and oxidation in response to changing environmental conditions, such as temperature and diet".


        # 221210: WARS looked very salt and pepper in the enrichment view in mcview, and due to its presumable function (note that https://www.proteinatlas.org/ENSG00000140105-WARS1 says "is induced by interferon"), decided to mark it as lateral.
        'WARS', # https://www.proteinatlas.org/ENSG00000140105-WARS1: "Tryptophanyl-tRNA synthetase 1" and "Aminoacyl-tRNA synthetases catalyze the aminoacylation of tRNA by their cognate amino acid. Because of their central role in linking amino acids with nucleotide triplets contained in tRNAs, aminoacyl-tRNA synthetases are thought to be among the first proteins that appeared in evolution. Two forms of tryptophanyl-tRNA synthetase exist, a cytoplasmic form, named WARS, and a mitochondrial form, named WARS2. Tryptophanyl-tRNA synthetase (WARS) catalyzes the aminoacylation of tRNA(trp) with tryptophan and is induced by interferon."

        
        # 'BCL2', # https://www.proteinatlas.org/ENSG00000171791-BCL2: "BCL2 apoptosis regulator" and "Suppresses apoptosis in a variety of cell systems including factor-dependent lymphohematopoietic and neural cells 1, 2. Regulates cell death by controlling the mitochondrial membrane permeability 3. Appears to function in a feedback loop system with caspases"

        
        # 'BCL2A1', # https://www.proteinatlas.org/ENSG00000140379-BCL2A1: "BCL2 related protein A1" and "Retards apoptosis induced by IL-3 deprivation. May function in the response of hemopoietic cells to external signals and in maintaining endothelial survival during infection (By similarity). Can inhibit apoptosis induced by serum starvation in the mammary epithelial cell line HC11 (By similarity)"

        # 221210: IFITM2 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'IFITM2', # https://www.proteinatlas.org/ENSG00000185201-IFITM2: "Interferon induced transmembrane protein 2" and "IFN-induced antiviral protein which inhibits the entry of viruses to the host cell cytoplasm, permitting endocytosis, but preventing subsequent viral fusion and release of viral contents into the cytosol 1, 2. Active against multiple viruses"

        'CD52', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, clustered somewhat nicely with IFITM2 and multiple HLA-D* genes. https://www.proteinatlas.org/ENSG00000169442-CD52.


        # 221210: IFIH1 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'IFIH1', # https://www.proteinatlas.org/ENSG00000115267-IFIH1: "Interferon induced with helicase C domain 1" and "Innate immune receptor which acts as a cytoplasmic sensor of viral nucleic acids and plays a major role in sensing viral infection and in the activation of a cascade of antiviral responses including the induction of type I interferons and pro-inflammatory cytokines"

        # 221210: IFIT5 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'IFIT5', # https://www.proteinatlas.org/ENSG00000152778-IFIT5: "Interferon induced protein with tetratricopeptide repeats 5" and "Interferon-induced RNA-binding protein involved in the human innate immune response. Has a broad and adaptable RNA structure recognition important for RNA recognition specificity in antiviral defense."

        # 221210: IFI35 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'IFI35', # https://www.proteinatlas.org/ENSG00000068079-IFI35: "Interferon induced protein 35" and "Acts as a signaling pathway regulator involved in innate immune system response 1, 2, 3. In response to interferon IFN-alpha, associates in a complex with signaling pathway regulator NMI to regulate immune response"

        'NMI', # https://www.proteinatlas.org/ENSG00000123609-NMI: "N-myc and STAT interactor" and "Acts as a signaling pathway regulator involved in innate immune system response" and "In response to interferon IFN-alpha, associates in a complex with signaling pathway regulator IFI35 to regulate immune response; the complex formation prevents proteasome-mediated degradation of IFI35 8, 9. In complex with IFI35, inhibits virus-triggered type I IFN-beta production when ubiquitinated by ubiquitin-protein ligase TRIM21 10. In complex with IFI35, negatively regulates nuclear factor NF-kappa-B signaling by inhibiting the nuclear translocation, activation and transcription of NF-kappa-B subunit p65/RELA, resulting in the inhibition of endothelial cell proliferation, migration and re-endothelialization of injured arteries 11. Negatively regulates virus-triggered type I interferon/IFN production by inducing proteosome-dependent degradation of IRF7, a transcriptional regulator of type I IFN, thereby interfering with cellular antiviral responses (By similarity). Beside its role as an intracellular signaling pathway regulator, also functions extracellularly as damage-associated molecular patterns (DAMPs) to promote inflammation, when actively released by macrophage to the extracellular space during cell injury or pathogen invasion"

        # 221210: IFI30 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'IFI30', # https://www.proteinatlas.org/ENSG00000216490-IFI30: "IFI30 lysosomal thiol reductase" and "Lysosomal thiol reductase that can reduce protein disulfide bonds." and https://www.genecards.org/cgi-bin/carddisp.pl?gene=IFI30: "Gamma-Interferon-Inducible Lysosomal Thiol Reductase"
        
        'LINC01857', # 230221: was highly correlated with IFI30, H2AFY, CD24. also, seems like it should be active only in B cells (in our data).

        # 221210: IFI30 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'IFITM3', # https://www.proteinatlas.org/ENSG00000142089-IFITM3: "Interferon induced transmembrane protein 3" and "IFN-induced antiviral protein which disrupts intracellular cholesterol homeostasis. Inhibits the entry of viruses to the host cell cytoplasm by preventing viral fusion with cholesterol depleted endosomes. May inactivate new enveloped viruses which buds out of the infected cell, by letting them go out with a cholesterol depleted membrane. Active against multiple viruses"

        'TNFSF10', # highly correlated with IFITM3 (0.9). https://www.proteinatlas.org/ENSG00000121858-TNFSF10: "TNF superfamily member 10"

        # 221210: IFI27 looked very salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'IFI27', # https://www.proteinatlas.org/ENSG00000165949-IFI27: "Interferon alpha inducible protein 27" and "Probable adapter protein involved in different biological processes 1, 2. Part of the signaling pathways that lead to apoptosis 3, 4, 5. Involved in type-I interferon-induced apoptosis characterized by a rapid and robust release of cytochrome C from the mitochondria and activation of BAX and caspases 2, 3, 6, 8 and 9"


        
        # 'IGKC', # looks batchy, but it seems actually more specific to donors than to batches, so not marking as lateral.

        
        # 'DLK1', # https://www.proteinatlas.org/ENSG00000185559-DLK1: "Delta like non-canonical Notch ligand 1" and "This gene encodes a transmembrane protein that contains multiple epidermal growth factor repeats that functions as a regulator of cell growth."

        'STK17B', # role in apoptosis? also, looked quite batchy in the kolmogorov smirnov tests when ignoring one donor of the experiment (though it was lower rather than higher, which is what you would expect for a specific experiment with more apoptosis). https://www.proteinatlas.org/ENSG00000081320-STK17B: "Serine/threonine kinase 17b" and "Phosphorylates myosin light chains (By similarity). Acts as a positive regulator of apoptosis."

        'PMAIP1', # https://www.proteinatlas.org/ENSG00000141682-PMAIP1: "Phorbol-12-myristate-13-acetate-induced protein 1" and "Promotes activation of caspases and apoptosis. Promotes mitochondrial membrane changes and efflux of apoptogenic proteins from the mitochondria. Contributes to p53/TP53-dependent apoptosis after radiation exposure."

        # 221210: XAF1 looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral.
        'XAF1', # https://www.proteinatlas.org/ENSG00000132530-XAF1: "XIAP associated factor 1" and "Seems to function as a negative regulator of members of the IAP (inhibitor of apoptosis protein) family. Inhibits anti-caspase activity of BIRC4. Induces cleavage and inactivation of BIRC4 independent of caspase activation. Mediates TNF-alpha-induced apoptosis and is involved in apoptosis in trophoblast cells. May inhibit BIRC4 indirectly by activating the mitochondrial apoptosis pathway. After translocation to mitochondria, promotes translocation of BAX to mitochondria and cytochrome c release from mitochondria."


        
        'HSP90AB1', # https://www.proteinatlas.org/ENSG00000096384-HSP90AB1: "Heat shock protein 90 alpha family class B member 1" and "Molecular chaperone that promotes the maturation, structural maintenance and proper regulation of specific target proteins involved for instance in cell cycle control and signal transduction." and "Biological process (UniProt)i Stress response"

        # decided to not mark as lateral among others because https://www.proteinatlas.org/ENSG00000118513-MYB says "This protein plays an essential role in the regulation of hematopoiesis."
        # 'MYB', # correlated with HSP90AB1, but didn't look batchy. https://www.proteinatlas.org/ENSG00000118513-MYB: "MYB proto-oncogene, transcription factor" and "This protein plays an essential role in the regulation of hematopoiesis."

        # 221210: HSP90AA1 is a highly expressed gene and a strong feature gene currently (IIUC). also, its presumable function seems relatively lateral. also looks a bit batchy.
        'HSP90AA1', # https://www.proteinatlas.org/ENSG00000080824-HSP90AA1: "Heat shock protein 90 alpha family class A member 1" and "Molecular chaperone that promotes the maturation, structural maintenance and proper regulation of specific target proteins involved for instance in cell cycle control and signal transduction."

        # listed in suspect_gene_names in https://github.com/tanaylab/metacells/blob/master/vignettes/Metacells_Vignette.ipynb
        'HSPA1A', # https://www.proteinatlas.org/ENSG00000204389-HSPA1A: "Heat shock protein family A (Hsp70) member 1A" and "Molecular chaperone implicated in a wide variety of cellular processes, including protection of the proteome from stress, folding and transport of newly synthesized polypeptides, activation of proteolysis of misfolded proteins and the formation and dissociation of protein complexes"

        'HSPE1', # https://www.proteinatlas.org/ENSG00000115541-HSPE1: "Heat shock protein family E (Hsp10) member 1" and "Co-chaperonin implicated in mitochondrial protein import and macromolecular assembly. Together with Hsp60, facilitates the correct folding of imported proteins. May also prevent misfolding and promote the refolding and proper assembly of unfolded polypeptides generated under stress conditions in the mitochondrial matrix"

        'HSP90B1', # marking as lateral mainly due to its presumable lateral function. also, looks (maybe) somewhat batchy in /dummy/dummy/dummy/raid/mds/metacell_analysis_only_ultima/intermediate_output/potential_batch_DEGs/downsampled_HSP90B1_UMI_cell_cumulative_hist_per_exp_and_donor_heatmap.png, and a bit salt and pepper in my current model. https://www.proteinatlas.org/ENSG00000166598-HSP90B1: "Heat shock protein 90 beta family member 1" and "Molecular chaperone that functions in the processing and transport of secreted proteins (By similarity)".

        'HSPA5', # marking as lateral mainly due to its presumable lateral function. also, looks somewhat batchy in /dummy/dummy/dummy/raid/mds/metacell_analysis_only_ultima/intermediate_output/potential_batch_DEGs/downsampled_HSPA5_UMI_cell_cumulative_hist_per_exp_and_donor_heatmap.png, and salt and pepper in my current model. https://www.proteinatlas.org/ENSG00000044574-HSPA5: "Heat shock protein family A (Hsp70) member 5" and "Endoplasmic reticulum chaperone that plays a key role in protein folding and quality control in the endoplasmic reticulum lumen 1, 2, 3, 4. Involved in the correct folding of proteins and degradation of misfolded proteins via its interaction with DNAJC10/ERdj5, probably to facilitate the release of DNAJC10/ERdj5 from its substrate (By similarity). Acts as a key repressor of the ERN1/IRE1-mediated unfolded protein response (UPR)" and "Biological process (UniProt)i Host-virus interaction"

        'HSPA8', # marking as lateral mainly due to its presumable lateral function. also, looks (maybe) somewhat batchy in /dummy/dummy/dummy/raid/mds/metacell_analysis_only_ultima/intermediate_output/potential_batch_DEGs/downsampled_HSPA8_UMI_cell_cumulative_hist_per_exp_and_donor_heatmap.png, and salt and pepper in my current model. https://www.proteinatlas.org/ENSG00000109971-HSPA8: "Heat shock protein family A (Hsp70) member 8" and "Molecular chaperone implicated in a wide variety of cellular processes, including protection of the proteome from stress, folding and transport of newly synthesized polypeptides, activation of proteolysis of misfolded proteins and the formation and dissociation of protein complexes. Plays a pivotal role in the protein quality control system, ensuring the correct folding of proteins, the re-folding of misfolded proteins and controlling the targeting of proteins for subsequent degradation 1, 2, 3, 4, 5. This is achieved through cycles of ATP binding, ATP hydrolysis and ADP release, mediated by co-chaperones"

        'HSPB1', # marking as lateral mainly due to its presumable lateral function. also, looks (maybe) somewhat batchy in /dummy/dummy/dummy/raid/mds/metacell_analysis_only_ultima/intermediate_output/potential_batch_DEGs/downsampled_HSPB1_UMI_cell_cumulative_hist_per_exp_and_donor_heatmap.png, and salt and pepper in my current model. https://www.proteinatlas.org/ENSG00000106211-HSPB1: "Heat shock protein family B (small) member 1" and "Small heat shock protein which functions as a molecular chaperone probably maintaining denatured proteins in a folding-competent state 1, 2. Plays a role in stress resistance and actin organization 3. Through its molecular chaperone activity may regulate numerous biological processes including the phosphorylation and the axonal transport of neurofilament proteins" and "This gene encodes a member of the small heat shock protein (HSP20) family of proteins. In response to environmental stress, the encoded protein translocates from the cytoplasm to the nucleus and functions as a molecular chaperone that promotes the correct folding of other proteins. This protein plays an important role in the differentiation of a wide variety of cell types. Expression of this gene is correlated with poor clinical outcome in multiple human cancers, and the encoded protein may promote cancer cell proliferation and metastasis, while protecting cancer cells from apoptosis. Mutations in this gene have been identified in human patients with Charcot-Marie-Tooth disease and distal hereditary motor neuropathy".

        'HSPH1', 

        'HSPA14', # https://www.proteinatlas.org/ENSG00000187522-HSPA14: "Heat shock protein family A (Hsp70) member 14" and "Component of the ribosome-associated complex (RAC), a complex involved in folding or maintaining nascent polypeptides in a folding-competent state. In the RAC complex, binds to the nascent polypeptide chain, while DNAJC2 stimulates its ATPase activity"

        'HSPD1', # marking as lateral mainly due to its presumable lateral function. also, looks somewhat batchy in /dummy/dummy/dummy/raid/mds/metacell_analysis_only_ultima/intermediate_output/potential_batch_DEGs/downsampled_HSPD1_UMI_cell_cumulative_hist_per_exp_and_donor_heatmap.png, and salt and pepper in my current model. https://www.proteinatlas.org/ENSG00000144381-HSPD1: "Heat shock protein family D (Hsp60) member 1" and "This gene encodes a member of the chaperonin family. The encoded mitochondrial protein may function as a signaling molecule in the innate immune system. This protein is essential for the folding and assembly of newly imported proteins in the mitochondria.".

        # DNAJ, aka HSP40, https://en.wikipedia.org/wiki/Chaperone_DnaJ
        'DNAJC11', 'DNAJC16', 'DNAJC8', 'DNAJC6', 'DNAJB4', 'DNAJC27', 'DNAJC27-AS1', 'DNAJC5G', 'DNAJC10', 'DNAJB2', 'DNAJB8', 'DNAJB8-AS1', 'DNAJC13', 'DNAJC19', 'DNAJB11', 'DNAJB14', 'DNAJC21', 'DNAJC18', 'DNAJC30', 'DNAJC2', 'DNAJB9', 'DNAJB6', 'DNAJC5B', 'DNAJA1', 'DNAJB5-DT', 'DNAJB5', 'DNAJC25', 'DNAJC25-GNG10', 'DNAJC24', 'DNAJC4', 'DNAJB13', 'DNAJC1', 'DNAJC12', 'DNAJB12', 'DNAJC9', 'DNAJC9-AS1', 'DNAJC22', 'DNAJC14', 'DNAJC15', 'DNAJC3-DT', 'DNAJC3', 'DNAJC17', 'DNAJA4', 'DNAJA3', 'DNAJA2', 'DNAJC7', 'DNAJC5', 'DNAJB1', 'DNAJB7', 'DNAJC28',

        'B2M', 

        'GBP5', # https://www.proteinatlas.org/ENSG00000154451-GBP5: "Guanylate binding protein 5" and "As an activator of NLRP3 inflammasome assembly, plays a role in innate immunity and inflammation. Promotes selective NLRP3 inflammasome assembly in response to microbial and soluble, but not crystalline, agents".

        'CAT', # https://www.proteinatlas.org/ENSG00000121691-CAT: "Catalase" and "Occurs in almost all aerobically respiring organisms and serves to protect cells from the toxic effects of hydrogen peroxide. Promotes growth of cells including T-cells, B-cells, myeloid leukemia cells, melanoma cells, mastocytoma cells and normal and transformed fibroblast cells"


        'EGR1', # https://www.proteinatlas.org/ENSG00000120738-EGR1: "Early growth response 1" and "Transcriptional regulator 1. Recognizes and binds to the DNA sequence 5'-GCG(T/G)GGGCG-3'(EGR-site) in the promoter region of target genes (By similarity). Binds double-stranded target DNA, irrespective of the cytosine methylation status 2, 3. Regulates the transcription of numerous target genes, and thereby plays an important role in regulating the response to growth factors, DNA damage, and ischemia. Plays a role in the regulation of cell survival, proliferation and cell death. Activates expression of p53/TP53 and TGFB1, and thereby helps prevent tumor formation. Required for normal progress through mitosis and normal proliferation of hepatocytes after partial hepatectomy. Mediates responses to ischemia and hypoxia;"

        'SAMHD1', # https://www.proteinatlas.org/ENSG00000101347-SAMHD1: "SAM and HD domain containing deoxynucleoside triphosphate triphosphohydrolase 1" and "Protein that acts both as a host restriction factor involved in defense response to virus and as a regulator of DNA end resection at stalled replication forks"

        'AEN', # https://www.proteinatlas.org/ENSG00000181026-AEN: "Apoptosis enhancing nuclease" and "Exonuclease with activity against single- and double-stranded DNA and RNA. Mediates p53-induced apoptosis. When induced by p53 following DNA damage, digests double-stranded DNA to form single-stranded DNA and amplifies DNA damage signals, leading to enhancement of apoptosis".

        'HNRNPA2B1', # https://www.proteinatlas.org/ENSG00000122566-HNRNPA2B1: "Heterogeneous nuclear ribonucleoprotein A2/B1" and "Heterogeneous nuclear ribonucleoprotein (hnRNP) that associates with nascent pre-mRNAs, packaging them into hnRNP particles. The hnRNP particle arrangement on nascent hnRNA is non-random and sequence-dependent and serves to condense and stabilize the transcripts and minimize tangling and knotting. Packaging plays a role in various processes such as transcription, pre-mRNA processing, RNA nuclear export, subcellular location, mRNA translation and stability of mature mRNAs" and "Also plays a role in the activation of the innate immune response 7. Mechanistically, senses the presence of viral DNA in the nucleus, homodimerizes and is demethylated by JMJD6 8. In turn, translocates to the cytoplasm where it activates the TBK1-IRF3 pathway, leading to interferon alpha/beta production".

        'PRKDC', # seems correlated with PCNA. https://www.proteinatlas.org/ENSG00000253729-PRKDC: "Protein kinase, DNA-activated, catalytic subunit" and "Serine/threonine-protein kinase that acts as a molecular sensor for DNA damage 1, 2, 3, 4. Involved in DNA non-homologous end joining (NHEJ) required for double-strand break (DSB) repair and V(D)J recombination"

        'GADD45A', # https://www.proteinatlas.org/ENSG00000116717-GADD45A: "Growth arrest and DNA damage inducible alpha" and "In T-cells, functions as a regulator of p38 MAPKs by inhibiting p88 phosphorylation and activity (By similarity). Might affect PCNA interaction with some CDK (cell division protein kinase) complexes; stimulates DNA excision repair in vitro and inhibits entry of cells into S phase" and "This gene is a member of a group of genes whose transcript levels are increased following stressful growth arrest conditions and treatment with DNA-damaging agents. The protein encoded by this gene responds to environmental stresses by mediating activation of the p38/JNK pathway via MTK1/MEKK4 kinase. The DNA damage-induced transcription of this gene is mediated by both p53-dependent and -independent mechanisms"

        'BTG3', # https://www.proteinatlas.org/ENSG00000154640-BTG3: "BTG anti-proliferation factor 3" and "Overexpression impairs serum-induced cell cycle progression from the G0/G1 to S phase". also, was upregulated in cells that waited for 24/36/48h.

        'CD69', # was upregulated in cells that waited for 24/36/48h. https://www.proteinatlas.org/ENSG00000110848-CD69: "Involved in lymphocyte proliferation and functions as a signal transmitting receptor in lymphocytes, natural killer (NK) cells, and platelets"

        'DDIT3', # was upregulated in cells that waited for 24/36/48h. https://www.proteinatlas.org/ENSG00000175197-DDIT3: "Multifunctional transcription factor in endoplasmic reticulum (ER) stress response 1, 2, 3. Plays an essential role in the response to a wide variety of cell stresses and induces cell cycle arrest and apoptosis in response to ER stress"
        'DDIT4', # was upregulated in cells that waited for 24/36/48h. https://www.proteinatlas.org/ENSG00000168209-DDIT4: "DNA damage inducible transcript 4" and "Regulates cell growth, proliferation and survival via inhibition of the activity of the mammalian target of rapamycin complex 1 (mTORC1). Inhibition of mTORC1 is mediated by a pathway that involves DDIT4/REDD1, AKT1, the TSC1-TSC2 complex and the GTPase RHEB. Plays an important role in responses to cellular energy levels and cellular stress, including responses to hypoxia and DNA damage. Regulates p53/TP53-mediated apoptosis in response to DNA damage via its effect on mTORC1 activity. Its role in the response to hypoxia depends on the cell type; it mediates mTORC1 inhibition in fibroblasts and thymocytes, but not in hepatocytes (By similarity). Required for mTORC1-mediated defense against viral protein synthesis and virus replication (By similarity)."

        'HIF1A', # was upregulated in cells that waited for 24/36/48h. https://www.proteinatlas.org/ENSG00000100644-HIF1A: "Hypoxia inducible factor 1 subunit alpha" and "Functions as a master transcriptional regulator of the adaptive response to hypoxia".
        'ATF3', # was upregulated in cells that waited for 24/36/48h. https://www.proteinatlas.org/ENSG00000162772-ATF3: "Activating transcription factor 3" and "This gene encodes a member of the mammalian activation transcription factor/cAMP responsive element-binding (CREB) protein family of transcription factors. This gene is induced by a variety of signals, including many of those encountered by cancer cells, and is involved in the complex process of cellular stress response".
        'G0S2', # was upregulated in cells that waited for 24/36/48h. https://www.proteinatlas.org/ENSG00000123689-G0S2: "G0/G1 switch 2" and "Promotes apoptosis by binding to BCL2, hence preventing the formation of protective BCL2-BAX heterodimers".
        
        'UCP2', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and somewhat low kruskal wallis pval, and its function sounds lateral: https://www.proteinatlas.org/ENSG00000175567-UCP2: "Uncoupling protein 2" and "Mitochondrial uncoupling proteins (UCP) are members of the larger family of mitochondrial anion carrier proteins (MACP). UCPs separate oxidative phosphorylation from ATP synthesis with energy dissipated as heat, also referred to as the mitochondrial proton leak. UCPs facilitate the transfer of anions from the inner to the outer mitochondrial membrane and the return transfer of protons from the outer to the inner mitochondrial membrane. They also reduce the mitochondrial membrane potential in mammalian cells",

        'STK4', # https://www.proteinatlas.org/ENSG00000101109-STK4: "Serine/threonine kinase 4" and "Stress-activated, pro-apoptotic kinase which, following caspase-cleavage, enters the nucleus and induces chromatin condensation followed by internucleosomal DNA fragmentation. Key component of the Hippo signaling pathway which plays a pivotal role in organ size control and tumor suppression by restricting proliferation and promoting apoptosis". when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval, and high in 24h/36h/48h samples.
        'LEPROT', # when looking at donor diff expr across MEBEMP-M cells, clustered nicely with STK4.
        'PDZD8', # when looking at donor diff expr across MEBEMP-M cells, clustered nicely with STK4.

        'TRIM8', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval, and high in 24h/36h/48h samples. https://www.proteinatlas.org/ENSG00000171206-TRIM8: "Tripartite motif containing 8" and "E3 ubiquitin-protein ligase that participates in multiple biological processes including cell survival, differentiation, apoptosis, and in particular, the innate immune response 1, 2. Participates in the activation of interferon-gamma signaling by promoting proteasomal degradation of the repressor SOCS1 3. Plays a positive role in the TNFalpha and IL-1beta signaling pathways".
        'RNF11', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval, and high in 24h/36h/48h samples. https://www.proteinatlas.org/ENSG00000123091-RNF11: "Ring finger protein 11" and "Essential component of a ubiquitin-editing protein complex, comprising also TNFAIP3, ITCH and TAX1BP1, that ensures the transient nature of inflammatory signaling pathways".
        'IER5', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval, and high in 24h/36h/48h samples. https://www.proteinatlas.org/ENSG00000162783-IER5: "Immediate early response 5" and "Mediates positive transcriptional regulation of several chaperone genes during the heat shock response in a HSF1-dependent manner" and "Involved in the regulation of cell proliferation and resistance to thermal stress 8, 9, 10. Involved in the cell cycle checkpoint and survival in response to ionizing radiation".
        'PIM3', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval, and high in 24h/36h/48h samples. https://www.proteinatlas.org/ENSG00000198355-PIM3: "Pim-3 proto-oncogene, serine/threonine kinase" and "Proto-oncogene with serine/threonine kinase activity that can prevent apoptosis, promote cell survival and protein translation. May contribute to tumorigenesis through: the delivery of survival signaling through phosphorylation of BAD which induces release of the anti-apoptotic protein Bcl-X(L), the regulation of cell cycle progression, protein synthesis and by regulation of MYC transcriptional activity".

        # 'SPATS2L', # when looking at donor diff expr across all cells, then calculating gene-gene correlations and clustering genes, was in a cluster of many IFN genes, including ISG15. # TODO: can we see it also only in MEBEMP-M, or some other state/trajectory?
        # 'ANXA2', # when looking at donor diff expr across all cells, then calculating gene-gene correlations and clustering genes, was in a cluster of many IFN genes, including ISG15. # TODO: can we see it also only in MEBEMP-M, or some other state/trajectory?

        # *[
        #     # the following are here mainly (or only) due to their relation to Dissimilar-high-IFN:
        #     # 'H3F3B', # was highly correlated with ISG15 in Dissimilar-high-IFN metacells. https://www.proteinatlas.org/ENSG00000132475-H3-3B: "H3.3 histone B" and "Variant histone H3 which replaces conventional H3 in a wide range of nucleosomes in active genes. Constitutes the predominant form of histone H3 in non-dividing cells and is incorporated into chromatin independently of DNA synthesis. Deposited at sites of nucleosomal displacement throughout transcribed genes, suggesting that it represents an epigenetic imprint of transcriptionally active chromatin."

        #     # 'NCOA7', # correlated with JUN and FOS in Dissimilar-high-IFN metacells. https://www.proteinatlas.org/ENSG00000111912-NCOA7: "Nuclear receptor coactivator 7" and "Enhances the transcriptional activities of several nuclear receptors".
        #     'XRN1', # correlated with ADAR, SAMD9, IFI35. was highly correlated with ISG15 in Dissimilar-high-IFN metacells. https://www.proteinatlas.org/ENSG00000114127-XRN1: "5'-3' exoribonuclease 1" and "Major 5'-3' exoribonuclease involved in mRNA decay. Required for the 5'-3'-processing of the G4 tetraplex-containing DNA and RNA substrates. The kinetic of hydrolysis is faster for G4 RNA tetraplex than for G4 DNA tetraplex and monomeric RNA tetraplex. Binds to RNA and DNA (By similarity). Plays a role in replication-dependent histone mRNA degradation."
        #     'GBP2', # separated well between Dissimilar-high-IFN and normal metacells. https://www.proteinatlas.org/ENSG00000162645-GBP2: "Guanylate binding protein 2" and "Hydrolyzes GTP to GMP in 2 consecutive cleavage reactions, but the major reaction product is GDP 1. Exhibits antiviral activity against influenza virus. Promotes oxidative killing and delivers antimicrobial peptides to autophagolysosomes, providing broad host protection against different pathogen classes (By similarity). Confers protection to the protozoan pathogen Toxoplasma gondii (By similarity)".

        #     'TYMP', # separated well between Dissimilar-high-IFN and normal metacells. https://www.proteinatlas.org/ENSG00000025708-TYMP: "Thymidine phosphorylase" and "May have a role in maintaining the integrity of the blood vessels. Has growth promoting activity on endothelial cells, angiogenic activity in vivo and chemotactic activity on endothelial cells in vitro". https://www.pnas.org/doi/10.1073/pnas.2105390118: "TYMP is an interferon-regulated gene in human"

        #     # the following were marked as lateral mainly because they separated Dissimilar-high-IFN metacells from other metacells, but were also diff expressed within Dissimilar-high-IFN, when comparing N208 cells to cells of others (controlling for metacell cluster). i am not sure we really want to mark them as lateral.
        #     # (more genes that (i think) were higher in N208 (compared in the same manner): ['RHEX', 'CD164', 'GBP4'])
        # https://www.stemcell.com/human-hematopoietic-stem-and-progenitor-cell-phenotyping-panels.html: Given table 1, marking cd38 as lateral doesn't seem like a good idea?
        #     'CD38', # separated well between Dissimilar-high-IFN and normal metacells. https://www.proteinatlas.org/ENSG00000004468-CD38
        #     'PPM1K', # separated well between Dissimilar-high-IFN and normal metacells, and was correlated with ISG15 in Dissimilar-high-IFN metacells. https://www.proteinatlas.org/ENSG00000163644-PPM1K: "Protein phosphatase, Mg2+/Mn2+ dependent 1K". https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2947462/: listed in table 1 ("Significantly Up-Regulated Human Hepatocyte Genes by IFN-α and IFN-γ at 6 and 18 h")
        #     'PDE7A', # # separated well between Dissimilar-high-IFN and normal metacells. https://www.proteinatlas.org/ENSG00000205268-PDE7A: "Phosphodiesterase 7A" and "Hydrolyzes the second messenger cAMP, which is a key regulator of many important physiological processes 1, 2, 3. May have a role in muscle signal transduction".
        #     'CHMP5', # the 4 genes it was most correlated with (considering only metacells of interest, i.e., non-disease and stem) were ISG15, IFIT3, OASL and IFIT2. https://www.proteinatlas.org/ENSG00000086065-CHMP5: "Charged multivesicular body protein 5" and "Probable peripherally associated component of the endosomal sorting required for transport complex III (ESCRT-III) which is involved in multivesicular bodies (MVBs) formation and sorting of endosomal cargo proteins into MVBs. MVBs contain intraluminal vesicles (ILVs) that are generated by invagination and scission from the limiting membrane of the endosome and mostly are delivered to lysosomes enabling degradation of membrane proteins, such as stimulated growth factor receptors, lysosomal enzymes and lipids".
        # ],
    ],

    # 'misc_batchy': [ # 240421: NOTE: commenting out for now, because lateral genes should belong to clear lateral programs etc (e.g., platelets, interferon).
        

    #     'EIF1AX', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, seemed somewhat batchy, and clearly lower in the batchy 07_02_22 experiments (N223, N224, N225, N226). also, its function sounds lateral: https://www.proteinatlas.org/ENSG00000173674-EIF1AX: "Eukaryotic translation initiation factor 1A X-linked" and "Seems to be required for maximal rate of protein biosynthesis. Enhances ribosome dissociation into subunits and stabilizes the binding of the initiator Met-tRNA(I) to 40 S ribosomal subunits".
    #     'DEK', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, seemed somewhat batchy, and clearly lower in the batchy 07_02_22 experiments (N223, N224, N225, N226).
        
    #     'NUCKS1', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, seemed batchy. also, its function sounds lateral: https://www.proteinatlas.org/ENSG00000069275-NUCKS1: "Nuclear casein kinase and cyclin dependent kinase substrate 1" and "Chromatin-associated protein involved in DNA repair by promoting homologous recombination (HR)".

    #     'HMGN2', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, clustered nicely with multiple genes that looked batchy.
    #     # 'FABP5', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, clustered nicely with multiple genes that looked batchy. but kruskal wallis pval was not low enough.
    #     # 'SNHG7', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, seemed somewhat batchy, and clustered nicely with VIM (which looked batchy). also, RNA gene. https://genome.ucsc.edu/cgi-bin/hgGene?hgg_gene=ENST00000414282.5&hgg_chrom=chr9&hgg_start=136721365&hgg_end=136728184&hgg_type=knownGene&db=hg38: "small nucleolar RNA host gene 7", but kruskal wallis pval was not low enough.
    #     'NARS', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, clustered nicely with multiple genes that looked batchy, and also looked batchy in a stripplot ordered by bleeding date mean log ratio.
    #     'RBMX', # when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, clustered nicely with multiple genes that looked batchy, and also looked batchy in a stripplot ordered by bleeding date mean log ratio.
    #     'NUCB2', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio.
    #     # 'CYTL1', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio. actually doesnt look that batchy in this plot... not very tight per bleeidng date... see the same plot for NUCB2 for comparison.
    #     'TFPI', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'SMIM24', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'PBXIP1', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'LINC02573', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval, and high in 24h/36h/48h samples.
    #     'MSI2', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'ANKRD12', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'SQSTM1', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval, and high in 24h/36h/48h samples.
    #     'MED30', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval, and high in 24h/36h/48h samples.
    #     'TSPO', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'MRPL16', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'PRDX1', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'MLLT3', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'CAPZA2', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'TFDP2', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'IGFBP7', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'EBPL', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'SNX3', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'COX5A', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'MAP1LC3B', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'GABARAPL2', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'UBE2B', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'RAC1', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'PEBP1', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'HNRNPD', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'PDIA3', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.

    #     # the following genes, when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, clustered nicely with multiple genes that looked batchy in a stripplot ordered by bleeding date mean log ratio and had low kruskal wallis pval:
    #     'PCBP1',
    #     'SUMO3',
    #     'NFKBIA',
    #     'SARAF',
    #     'HNRNPA0',
    #     'UBE2D2',
    #     'UQCRFS1',
    #     'HEBP2',
    #     'YBX3',
    #     'LAPTM4B',
    #     'NFE2',
    #     'CCNG1',
    #     'EIF2S3',
    #     'FKBP1A',
    #     'MDH2',
    #     'HMGA1',
    #     'PHB2',
    #     'HEBP1',
    #     'SPINT2',
    #     'MARCKSL1',
    #     'PPDPF',
    #     'SDCBP',
    #     'TUBA1A',
    #     'SNRPN',
    #     'LDHA',
    #     'STRAP',
    #     'CNBP',
    #     'RHOA',
    #     'GLUL',
    #     'RHOG',
    #     'ARF1',
    #     'PSMD7',
    #     'ATF4',
    #     'ENSA',
    #     'UBC',
    #     'DDAH2',
    #     'LAPTM5',
    #     'CIRBP',
    #     'CLIC1',
    #     'CASP4',
    #     'BRK1',
    #     'MYL12A',
    #     'NPC2',
    #     'ATP6V1G1',
    #     'UBB',
    #     'LAPTM4A',
    #     'TMEM59',
    #     'RNF7',
    #     'PSMB1',
    #     'SAP18',
    #     'CCT8',
    #     'SSB',
    #     'HSD17B11',
    #     'RSL1D1',
    #     'HIGD2A',
    #     'ATP6V1F',
    #     'SSR2',
    #     'NDUFV2',
    #     'PARK7',
    #     'SRI',
    #     'SLC25A5',
    #     'APEX1',
    #     'GDI2',
    #     'ERP29',
    #     'PSMG2',
    #     'NDUFB5',
    #     'OSTC',
    #     'IGBP1',
    #     'UBXN1',
    #     'PRDX6',
    #     'EIF3H',
    #     'EIF3M',
    #     'RTRAF',
    #     'COX7A2L',
    #     'ATP5F1C',
    #     'PSMA1',
    #     'CSNK2B',
    #     'ILF2',
    #     'TXNL1',
    #     'ATP5PB',
    #     'CCT4',
    #     'OLA1',
    #     'VDAC2',
    #     'HADHB',
    #     'UQCRC2',
    #     'ATP5F1B',
    #     'EIF3I',
    #     'RAC2',
    #     'UBE2E3',
    #     'SYPL1',
    #     'TMBIM6',
    #     'ZNF22',

    #     'PRDX5', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'HMGN3', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'C11orf1', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
    #     'RWDD1', # when looking at donor diff expr across MEBEMP-M cells, looked batchy in a stripplot ordered by bleeding date mean log ratio, and low kruskal wallis pval.
        
    #     # the following genes, when looking at donor diff expr across MEBEMP-M cells and clustering genes accordingly, clustered nicely with multiple genes that looked batchy in a stripplot ordered by bleeding date mean log ratio and had low kruskal wallis pval:
    #     'SUB1',
    #     'HSBP1',
    #     'NDUFC2',
    #     'SUMO2',
    #     'ESD',
    #     'SPCS1',
    #     'GSTO1',
    #     'MPC2',
    #     'PDCD5',
    #     'NDUFB6',
    #     'SNRPD1',
    #     'ATP5MC3',
    #     'NDUFA12',
    #     'ATP5IF1',
    #     'MINOS1',
    #     'SNRPG',
    #     'RBX1',
    #     'COX6C',
    #     'COX7A2',
    #     'ATP5PF',
    #     'NDUFA4',
    #     'SEC61B',
    #     'KRT10',
    #     'CHCHD2',
    #     'GNG5',
    #     'RBM3',
    # ],
    # 'batchy_by_kruskal': [ # NOTE: 240328: commented out after seeing GATA2 in UPREGULATED_IN_ANY_STATE_SET_IN_DELAYED_GENES.
    #     *BATCHY_IN_ANY_STATE_SET_BY_KRUSKAL_GENES,
    # ],
    # 'upregulated_in_samples_that_waited': [ # NOTE: 240328: commented out after seeing GATA2 in UPREGULATED_IN_ANY_STATE_SET_IN_DELAYED_GENES.
    #     *UPREGULATED_IN_ANY_STATE_SET_IN_DELAYED_GENES,
    # ],
    'misc_noisy': [ # i guess this mostly makes sense if we consider ambient noise? as there would always be "highest expressed" genes...
        'IGLC2', 

        'PRSS2', # 240424: see Dropbox\Tanay_MDS_etc\240403_donor_diff_exp\clp\prss2_sig\prss2_corr_with_rest_of_prss2_sig, specifically prss2_across_nimrod_atlas_c_clp_stratified_by_ds_prss2_sig_excluding_prss2.png. seems like PRSS2 is correlated with the rest of prss2_sig, yet it's variance seems high, which is what led to the 230623 observation, i guess. but why didn't we see on 230623 something similar for CLC (for example) in BEMPs? maybe because in BEMPs we have multiple such highly variable genes, while PRSS2 is the only such gene for CLP?? 230623: some CLP-M MCs have a 16-fold higher PRSS2 expr than others, without anything else changing much. so maybe it is some crazy expression of the T cell receptor locus? anyway, seems like it shouldn't be a feature gene due to this "noisiness". (maybe not a sufficient reason to mark as lateral, but note that PRSS2 is non-coding. And in T cell receptor beta locus (https://genome.ucsc.edu/cgi-bin/hgGene?hgg_gene=ENST00000539842.6&hgg_chrom=chr7&hgg_start=142770969&hgg_end=142774560&hgg_type=knownGene&db=hg38).)
        # 231119: NOTE: actually i don't want to mark these other T cell receptor genes as lateral, as i didn't notice noisiness in them, and they might actually tell me about whether the chromatin there is open, maybe???
        # 'TRBVB', # correlated with PRSS2 after marking PRSS2 noisy. https://www.genecards.org/cgi-bin/carddisp.pl?gene=TRBVB: "T Cell Receptor Beta Variable B (Pseudogene)"

        # 'MS4A2', # 231121: it seems to sometimes have a 20-fold higher expression in BEMP (or post-BEMP), yet it is correlated somewhat nicely with MS4A3 (in BEMP and post-BEMP), so don't mark as noisy...

        'NEAT1', # 221210: on the first metacell run, NEAT1 was the gene with highest 'feature_gene'. given it is not that reliable (due to its nucleus localization), I decided it is better to mark it as lateral.

        'S100A10', # S100A10 has a weird jump in MEBEMP-M and MPP. maybe bursty? anyway, it is quite highly correlated with ANXA2 (in each cell state), (somewhat) except for the jump, and so it seems to make sense to mark it as noisy, as we won't lose the information thanks to ANXA2. (also nicely correlated with S100A11.)
        
        'AC011139.1', # AC011139.1 (ENSG00000283458) has a weird jump in MEBEMP-M, MPP and HSC. maybe bursty? it is probably non-coding RNA (http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000283458;r=18:39314224-39328009: gene type: LncRNA). anyway, it is correlated with PHLDB2 (in each cell state), (somewhat) except for the jump, and so it seems to make sense to mark it as noisy, as we won't lose the information thanks to PHLDB2.

        'EREG', # was upregulated in cells that waited for 24/36/48h. there were also two MEBEMP-M and MEBEMP-L metacells with a much higher EREG expression than others, and in MEBEMP-M and MEBEMP-L EREG is a bit lonely. also, it was identified by MCView QC as a gene that often has high inner-fold in MCs...
        
        'HBD', # pretty much correlated with nothing. specifically see 

        'LTB', # pretty much correlated with nothing.

        # 221210: H1FX looked quite salt and pepper in the enrichment view in mcview, and due to its presumable function, decided to mark it as lateral. https://www.proteinatlas.org/ENSG00000184897-H1-10 says "This gene encodes a replication-independent histone". but i still think it is lateral - it is about maintenance of histones, IIUC, which i don't think is specific to a cell type (https://www.sciencedirect.com/science/article/pii/S2211124716313304?via%3Dihub: "Replication-independent histone variants control transcriptional plasticity in postmitotic cells" and "While canonical histones are highly expressed exclusively during S phase, variant counterparts are expressed throughout the cell cycle and their incorporation occurs in a replication-independent manner (Weber and Henikoff, 2014). Exciting new insights have uncovered the relevance of histone variants as regulatory elements of local chromatin structure and dynamics").
        # https://www.proteinatlas.org/ENSG00000184897-H1-10: "H1.10 linker histone" and "H1-10 (H1FX, H1X, MGC15959, MGC8350)" and "Histones H1 are necessary for the condensation of nucleosome chains into higher-order structures" and "Nucleosomes consist of approximately 146 bp of DNA wrapped around a histone octamer composed of pairs of each of the four core histones (H2A, H2B, H3, and H4). The chromatin fiber is further compacted through the interaction of a linker histone, H1, with the DNA between the nucleosomes to form higher order chromatin structures. This gene encodes a replication-independent histone that is a member of the histone H1 family."
        'H1FX', # 230716: H1FX had a jump in each of MPP, MEBEMP-M and HSC. maybe bursty? anyway, the jump was pretty lonely in each of these states. also, it is a histone gene, so feels relatively safe to mark as noisy. also, it was identified by MCView QC as a gene that often has high inner-fold in MCs...
        'H1F0', # was upregulated in cells that waited for 24/36/48h. also, it was identified by MCView QC as a gene that often has high inner-fold in MCs. also, a histone.

        *NOISY_6P_HIST_GENE_NAMES,

        # TODO: maybe make them not noisy if we use MCNoise, but i guess why not just keep marking them as noisy. they shouldn't make all the difference on their on...
        # the following seem really noisy. my guess is ambient noise and also droplets with variable amounts of RBCs/reticulocytes/etc, and then the metacell algorithm cluster according to this variable amount of ambient noise/RBCs/reticuloctyes in the droplet, which is obviously wrong. see why_forbid_hb_if_no_mcnoise.png. (on 230628, it seemed like an earlier model, that was created without skipping MCNoise, did not have this issue...)
        'HBB', # also was identified as bursty lonely on 240123 for a c_ad of all experiments with MPN patients.
        'HBA1',
        'HBA2',
        # similarly for monocytes/neutrophils/etc. see why_forbid_s100a8_s100a9_if_no_mcnoise.png.
        'S100A8',
        'S100A9',
        # i don't know why i don't see something similar for LYZ, even though it is almost as strong as S100A8. maybe it isn't as strong in neutrophils/monocytes/whatever.
    ],

    'highly_polymorphic_genes': [ # i.e., containing many common SNPs.
        # immunoglobulin variable region
        'IGKV4-1', 'IGKV5-2', 'IGKV7-3', 'IGKV1-5', 'IGKV1-6', 'IGKV3-7', 'IGKV1-8', 'IGKV1-9', 'IGKV2-10', 'IGKV3-11', 'IGKV1-12', 'IGKV3-15', 'IGKV1-16', 'IGKV1-17', 'IGKV2-18', 'IGKV3-20', 'IGKV6-21', 'IGKV2-24', 'IGKV2-26', 'IGKV1-27', 'IGKV2-29', 'IGKV2-30', 'IGKV1-37', 'IGKV1-39', 'IGKV1D-39', 'IGKV1D-37', 'IGKV1D-33', 'IGKV2D-30', 'IGKV2D-29', 'IGKV2D-28', 'IGKV2D-24', 'IGKV6D-21', 'IGKV3D-20', 'IGKV1D-17', 'IGKV1D-16', 'IGKV3D-15', 'IGKV1D-13', 'IGKV3D-11', 'IGKV1D-43', 'IGKV1D-8', 'IGKV2OR2-1', 'IGKV1OR2-3', 'IGKV1OR2-108', 'IGKV2OR22-3', 'IGKV2-4', 'IGKV3-34', 'IGKV1OR22-1', 'IGKV6D-41',
        'IGHV1-12', 'IGHV1-14', 'IGHV1-17', 'IGHV1-18', 'IGHV1-2', 'IGHV1-24', 'IGHV1-3', 'IGHV1-46', 'IGHV1-58', 'IGHV1-67', 'IGHV1-68', 'IGHV1-69', 'IGHV1-69-2', 'IGHV1-69D', 'IGHV1OR15-2', 'IGHV1OR15-3', 'IGHV1OR15-4', 'IGHV1OR15-6', 'IGHV1OR16-4', 'IGHV2-26', 'IGHV2-5', 'IGHV2-70', 'IGHV2-70D', 'IGHV2OR16-5', 'IGHV3-11', 'IGHV3-13', 'IGHV3-15', 'IGHV3-19', 'IGHV3-20', 'IGHV3-21', 'IGHV3-22', 'IGHV3-23', 'IGHV3-29', 'IGHV3-30', 'IGHV3-30-2', 'IGHV3-32', 'IGHV3-33', 'IGHV3-38', 'IGHV3-41', 'IGHV3-42', 'IGHV3-43', 'IGHV3-48', 'IGHV3-49', 'IGHV3-52', 'IGHV3-53', 'IGHV3-57', 'IGHV3-62', 'IGHV3-64', 'IGHV3-64D', 'IGHV3-65', 'IGHV3-66', 'IGHV3-69-1', 'IGHV3-7', 'IGHV3-71', 'IGHV3-72', 'IGHV3-73', 'IGHV3-74', 'IGHV3-75', 'IGHV3-76', 'IGHV3-79', 'IGHV3OR16-6', 'IGHV3OR16-7', 'IGHV3OR16-9', 'IGHV4-28', 'IGHV4-31', 'IGHV4-34', 'IGHV4-39', 'IGHV4-4', 'IGHV4-59', 'IGHV4-61', 'IGHV4-80', 'IGHV5-10-1', 'IGHV5-51', 'IGHV5-78', 'IGHV6-1', 'IGHV7-27', 'IGHV7-4-1', 'IGHV7-40', 'IGHVII-22-1', 'IGHVII-26-2', 'IGHVII-28-1', 'IGHVII-30-1', 'IGHVII-30-21', 'IGHVII-33-1', 'IGHVII-43-1', 'IGHVII-44-2', 'IGHVII-60-1', 'IGHVII-74-1', 'IGHVIII-2-1', 'IGHVIII-38-1', 'IGHVIII-44', 'IGHVIII-76-1', 'IGHVIV-44-1',
        'IGHV3-33-2', 'IGHV3-6', 'IGHV3-37', 'IGHV3OR16-8', 'IGHV1-45', 'IGHVII-78-1', 'IGHV3-63', 'IGHV3-35', 'IGHV3-25', 'IGHV4-55', 
        'IGHJ2P', 'IGHJ4', 'IGHJ6',
        'IGLJ3',
        'IGLV3-6', 'IGLV1-36', 'IGLV1-40', 'IGLV1-41', 'IGLV1-44', 'IGLV1-47', 'IGLV1-50', 'IGLV1-51', 'IGLV10-54', 'IGLV2-11', 'IGLV2-14', 'IGLV2-18', 'IGLV2-23', 'IGLV2-5', 'IGLV2-8', 'IGLV3-1', 'IGLV3-10', 'IGLV3-12', 'IGLV3-16', 'IGLV3-17', 'IGLV3-19', 'IGLV3-21', 'IGLV3-25', 'IGLV3-26', 'IGLV3-27', 'IGLV3-9', 'IGLV4-3', 'IGLV4-60', 'IGLV4-69', 'IGLV5-37', 'IGLV5-45', 'IGLV5-48', 'IGLV5-52', 'IGLV6-57', 'IGLV7-35', 'IGLV7-43', 'IGLV7-46', 'IGLV8-61', 'IGLV8OR8-1', 'IGLV9-49', 'IGLVI-56', 'IGLVI-70', 'IGLVIV-59', 'IGLVIVOR22-1', 'IGLVV-58', 'IGLV11-55', 'IGLV10-67', 

        # unordered - no time at the moment.
        'IGKV1-35', 'IGHV3OR16-10', 'IGLV3-22', 'IGKV2D-26', 'IGKV1D-12', 'IGKV2-28', 'IGKV2D-18', 'IGHV3-60', 'IGHVIII-82', 'IGKV1OR-3', 'IGKV1OR1-1', 'IGHV1OR21-1', 'IGHV7-81', 'IGKV2OR2-2', 'IGKV1D-35', 'IGHVII-65-1', 'IGKV2-23', 'IGLV3-2', 'IGHV1OR16-1', 'IGKV1-13', 'IGHV3OR15-7', 'IGHV3OR16-13', 'IGLV2-34', 'IGKV1D-27', 'IGHV3-16',

        # IG* - https://www.frontiersin.org/articles/10.3389/fimmu.2020.02016/full: "in spite of the “constant” naming convention, genes of the immunoglobulin heavy-chain constant (IGHC) and light-chain constant (IGLC) loci are polymorphic, although to a far lesser extent than the immunoglobulin heavy-chain variable locus"
        # https://onlinelibrary.wiley.com/doi/10.1111/j.1399-0039.2006.00729.x: "The kappa immunoglobulin light chain is polymorphic and restricted to B cells; hence, it has the potential to act as an mHA and be a target for the GVL phenomenon in patients with B-cell malignancies when the donor and recipient have different kappa alleles"
        'IGKC',
        # it is easy to identify B cells, so just mark as lateral all immunoglobulin genes.
        'IGHA2', 'IGHA1',
        'IGHD',
        'IGHE', 'IGHEP1', 'IGHEP2',
        'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHGP',
        'IGHM',
        'IGLC2', 'IGLC3', 'IGLC4', 'IGLC5', 'IGLC6', 'IGLC7',

        # according to https://en.wikipedia.org/wiki/T-cell_receptor, there are T cell receptor alpha, beta, gamma and delta loci.
        # i assume that they are also quite polymorphic, similarly to immunoglobulin genes. (i guess this paper supports this: T cell receptor beta germline variability is revealed by inference from repertoire data (https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-01008-4, 2022))
        'TRAV1-1', 'TRAV1-2', 'TRAV2', 'TRAV3', 'TRAV4', 'TRAV5', 'TRAV6', 'TRAV8-1', 'TRAV10', 'TRAV11', 'TRAV12-1', 'TRAV8-2', 'TRAV8-3', 'TRAV13-1', 'TRAV12-2', 'TRAV8-4', 'TRAV8-5', 'TRAV13-2', 'TRAV14DV4', 'TRAV9-2', 'TRAV12-3', 'TRAV8-6', 'TRAV16', 'TRAV17', 'TRAV18', 'TRAV19', 'TRAV20', 'TRAV21', 'TRAV22', 'TRAV23DV6', 'TRAV24', 'TRAV25', 'TRAV26-1', 'TRAV8-7', 'TRAV27', 'TRAV28', 'TRAV29DV5', 'TRAV30', 'TRAV31', 'TRAV32', 'TRAV33', 'TRAV26-2', 'TRAV34', 'TRAV35', 'TRAV36DV7', 'TRAV37', 'TRAV38-1', 'TRAV38-2DV8', 'TRAV39', 'TRAV40', 'TRAV41',
        'TRAC',
        'TRBV1', 'TRBV2', 'TRBV3-1', 'TRBV4-1', 'TRBV5-1', 'TRBV6-1', 'TRBV7-1', 'TRBV4-2', 'TRBV6-2', 'TRBV7-2', 'TRBV8-1', 'TRBV5-2', 'TRBV6-4', 'TRBV7-3', 'TRBV8-2', 'TRBV5-3', 'TRBV9', 'TRBV10-1', 'TRBV11-1', 'TRBV12-1', 'TRBV10-2', 'TRBV11-2', 'TRBV12-2', 'TRBV6-5', 'TRBV7-4', 'TRBV5-4', 'TRBV6-6', 'TRBV7-5', 'TRBV5-5', 'TRBV6-7', 'TRBV7-6', 'TRBV5-6', 'TRBV6-8', 'TRBV7-7', 'TRBV5-7', 'TRBV7-9', 'TRBV13', 'TRBV10-3', 'TRBV11-3', 'TRBV12-3', 'TRBV12-4', 'TRBV12-5', 'TRBV14', 'TRBV15', 'TRBV16', 'TRBV17', 'TRBV18', 'TRBV19', 'TRBV20-1', 'TRBV21-1', 'TRBV22-1', 'TRBV23-1', 'TRBV24-1', 'TRBV25-1', 'TRBVA', 'TRBVB', 'TRBV26', 'TRBV27', 'TRBV28', 'TRBV29-1', 'TRBV30', 'TRBV20OR9-2', 'TRBV21OR9-2', 'TRBV23OR9-2', 'TRBV26OR9-2', 'TRBV29OR9-2',
        'TRBC1', 'TRBC2',
        'TRGV11', 'TRGVB', 'TRGV10', 'TRGV9', 'TRGVA', 'TRGV8', 'TRGV7', 'TRGV6', 'TRGV5P', 'TRGV5', 'TRGV4', 'TRGV3', 'TRGV2', 'TRGV1',
        'TRGC2', 'TRGC1',
        'TRDV1', 'TRDV2', 'TRDJ1', 'TRDV3',
        'TRDC', 



        # https://pubmed.ncbi.nlm.nih.gov/11486048/: "The human leukocyte antigens (HLA) encoded by genes within the major histocompatibility complex display an impressive degree of polymorphism."
        # NOTE: also, while working on CNAs, it seemed that HLA genes are upregulated in most cells of demux_28_11_21_1 (i.e., batchy).
        'HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B', 'HLA-DRA', 'HLA-DRB5', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DQB1-AS1', 'HLA-DQA2', 'HLA-DQB2', 'HLA-DOB', 'HLA-DMB', 'HLA-DMA', 'HLA-DOA', 'HLA-DPA1', 'HLA-DPB1',

        # https://humgenomics.biomedcentral.com/articles/10.1186/s40246-018-0175-1: "Amongst > 20,000 genes in the human genome, beta haemoglobin (HBB) gene is the most polymorphic gene, containing approximately 176 SNVs per kilobase (kb) with the highest density of SNVs within its coding region (Fig. 3a, red) (570 SNVs/kb). Several other haemoglobin genes (in green boxes) are also amongst the most polymorphic genes in the human genome with the majority of their SNVs residing within coding exons (red). Other highly polymorphic genes include the MHC family of genes (blue box) with most of their SNVs residing within introns (Fig. 3a, green) as well as the olfactory receptor (OR) gene family (orange box) where all the SNVs are also found within the coding region (Fig 3a, red)."
        # TODO: mark HBB as lateral??

        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8306414/: "The ABO gene displays a substantial spectrum of single nucleotide polymorphisms (SNPs): over 100 unique alleles have been curated for the ABO gene in the Leiden Open Variation Database [4]." and "To investigate if COVID-19 severity was associated with specific SNPs located on the ABO locus, all the 240 SNPs located in the ABO candidate gene were extracted from a GWAS and analyzed in this study."
        'ABO',

        # TODO: this could be easily checked though - just look at cellsnp-lite output, right? (combined with BED for gene positions)

    ],

    'platelet_related': [
        *PROBABLY_PLATELET_SPECIFIC_GENES,
        *STRONGEST_PLATELET_GENES,

        # rationale (my reasoning):
        # we mark as lateral genes with strongest expression in platelets, assuming those will form most of the "ambient" noise due to platelets. i guess.

        # megakaryocyte and platelet stuff
        
        'PPBP', # https://www.proteinatlas.org/ENSG00000163736-PPBP: "Pro-platelet basic protein" and "LA-PF4 stimulates DNA synthesis, mitosis, glycolysis, intracellular cAMP accumulation, prostaglandin E2 secretion, and synthesis of hyaluronic acid and sulfated glycosaminoglycan. It also stimulates the formation and secretion of plasminogen activator by human synovial cells." and "in vitro released from activated platelet alpha-granules" and "This growth factor is a potent chemoattractant and activator of neutrophils." (max PF4-PPBP-high metacell expression on 230105: -4)
        
        'PF4', # https://www.proteinatlas.org/ENSG00000163737-PF4: "Platelet factor 4" and "Released during platelet aggregation." and "chemokine is released from the alpha granules of activated platelets" (max PF4-PPBP-high metacell expression on 230105: -5.4)
        'RGS18', # https://www.proteinatlas.org/ENSG00000150681-RGS18: "Regulator of G protein signaling 18" and "Single cell type expression clusteri Platelets - Platelet activation (mainly)" (mainly)" (max PF4-PPBP-high metacell expression on 230105: -5.7)
        'GNG11', # highly correlated with PF4. https://www.proteinatlas.org/ENSG00000127920-GNG11: "G protein subunit gamma 11" and "Single cell type expression clusteri Platelets - Platelet activation (mainly)" (max PF4-PPBP-high metacell expression on 230105: -6.1)
        'TUBB1', # https://www.proteinatlas.org/ENSG00000101162-TUBB1: "Tubulin beta 1 class VI" and "Tubulin is the major constituent of microtubules. It binds two moles of GTP, one at an exchangeable site on the beta chain and one at a non-exchangeable site on the alpha chain (By similarity)." and "Single cell type expression clusteri Platelets - Platelet activation (mainly)" (max PF4-PPBP-high metacell expression on 230105: -6.3)
        'RAP1B', # https://www.proteinatlas.org/ENSG00000127314-RAP1B: "RAP1B, member of RAS oncogene family" and "GTP-binding protein that possesses intrinsic GTPase activity." and "Single cell type expression clusteri Platelets - Platelet activation (mainly)" (max PF4-PPBP-high metacell expression on 230105: -6.3)
        'CAVIN2', # https://www.proteinatlas.org/ENSG00000168497-CAVIN2: "Caveolae associated protein 2" and "Single cell type expression clusteri Platelets - Platelet activation (mainly)" (max PF4-PPBP-high metacell expression on 230105: -6.7)
        'RGS10', # https://www.proteinatlas.org/ENSG00000148908-RGS10: "Regulator of G protein signaling 10" and "Single cell type expression clusteri Platelets - Platelet activation (mainly)" (mainly)" (max PF4-PPBP-high metacell expression on 230105: -7.3)
        'ACRBP', # highly correlated with GP9. https://www.proteinatlas.org/ENSG00000111644-ACRBP: "Acrosin binding protein" and "Single cell type expression clusteri Platelets - Platelet activation (mainly)". (max PF4-PPBP-high metacell expression on 230105: -7.6)
        'NRGN', # https://www.proteinatlas.org/ENSG00000154146-NRGN: "Neurogranin" and "Single cell type expression clusteri Platelets - Platelet activation (mainly)" (max PF4-PPBP-high metacell expression on 230105: -7.7)
        'TUBA4A', # https://www.proteinatlas.org/ENSG00000127824-TUBA4A: "Tubulin alpha 4a" and "Single cell type expression clusteri Platelets - Platelet activation (mainly)" (max PF4-PPBP-high metacell expression on 230105: -7.9)
        'PTCRA', # very highly correlated with ACRBP. https://www.proteinatlas.org/ENSG00000171611-PTCRA: "Pre T cell antigen receptor alpha" and "Single cell type expression clusteri	Platelets - Platelet activation (mainly)" (max PF4-PPBP-high metacell expression on 230105: -7.9)
        'TMEM40', # highly correlated with GP9. https://www.proteinatlas.org/ENSG00000088726-TMEM40: "Transmembrane protein 40" and "Single cell type expression clusteri Platelets - Platelet activation (mainly)". (max PF4-PPBP-high metacell expression on 230105: -8)
        'GP9', # highly correlated with PF4. https://www.proteinatlas.org/ENSG00000169704-GP9: "Glycoprotein IX platelet" and "This gene encodes a small membrane glycoprotein found on the surface of human platelets." (max PF4-PPBP-high metacell expression on 230105: -8.2)
        'MMD', # very highly correlated with ACRBP, PTCRA. https://www.proteinatlas.org/ENSG00000108960-MMD: "Monocyte to macrophage differentiation associated" and "Single cell type expression clusteri Platelets - Platelet activation (mainly)" (max PF4-PPBP-high metacell expression on 230105: -8.2)
        'TREML1', # https://www.proteinatlas.org/ENSG00000161911-TREML1: "Triggering receptor expressed on myeloid cells like 1" and "Cytoplasmic expression in megakaryocytes and subtypes of cells in the spleen with positivity in platelets in several tissues" and "Single cell type expression clusteri Platelets - Platelet activation (mainly)" and "This gene encodes a member of the triggering receptor expressed on myeloid cells-like (TREM) family. The encoded protein is a type 1 single Ig domain orphan receptor localized to the alpha-granule membranes of platelets." (max PF4-PPBP-high metacell expression on 230105: -8.6)
        'SH3BGRL3', # in metacells with very high PF4 and/or PPBP, SH3BGRL3 is also very high (in some metacells even more than -7 (log2)). https://www.proteinatlas.org/ENSG00000142669-SH3BGRL3: "SH3 domain binding glutamate rich protein like 3"
        'H3F3A', # in metacells with very high PF4 and/or PPBP, H3F3A is also very high (in some metacells even more than -7 (log2)). https://www.proteinatlas.org/ENSG00000163041-H3-3A: "H3.3 histone A" and "Variant histone H3 which replaces conventional H3 in a wide range of nucleosomes in active genes. Constitutes the predominant form of histone H3 in non-dividing cells and is incorporated into chromatin independently of DNA synthesis. Deposited at sites of nucleosomal displacement throughout transcribed genes, suggesting that it represents an epigenetic imprint of transcriptionally active chromatin."
        'GP1BA', # correlated with PF4. https://www.proteinatlas.org/ENSG00000185245-GP1BA: "Glycoprotein Ib platelet subunit alpha" and "GP-Ib, a surface membrane protein of platelets, participates in the formation of platelet plugs by binding to the A1 domain of vWF, which is already bound to the subendothelium."
        'TUBA8', # highly correlated with PF4. https://www.proteinatlas.org/ENSG00000183785-TUBA8: "Tubulin alpha 8"


        # more platelet/MK (at least associated) genes, identified mainly by correlation with other platelet/MK genes
        'GP6',

        'CLEC1B',
        'LGALSL',
        'CMTM5',
        'CLDN5',

        'CTTN',
        'PDE5A',
        'ESAM',

        'PDGFA',

        'LGALS12',
        'AQP10',
        'GNAZ',
        'TRAPPC3L',
        'HGD',
        'DMTN',
        'DENND2C',

        'SLFN14', # https://www.jci.org/articles/view/80347: "In humans and mice, SLFN14 is located in a SLFN cluster with other schlafen paralogs (9). Members of the SLFN family are highly conserved among mammalian species. SLFN family proteins contain a unique motif of unknown function, the “SLFN box,” and an AAA domain. The AAA+ domain consists of a P-loop NTPase implicated in ATP/GTP binding and hydrolysis (10). The SLFN family members are divided into 3 groups. SLFN5, SLFN8, SLFN9, SLFN10, and SLFN14 all belong to group 3, although SLFN14 is unique in containing a putative nuclear localization RKRRR motif in its C-terminus extension (10). The SLFN family of proteins have been suggested to be critical for a variety of processes, including cell-cycle regulation, proliferation, and differentiation (10–14). Recently, data have been published suggesting an important function for SLFN14 as an endoribonuclease, regulating rRNA and ribosome-associated mRNA cleavage and translational control in rabbit reticulocytes (15)."

        # the following genes were estimated to be at least 8-fold enriched in platelets, and with expression of at least -10 in platelets.
        'GRAP2', # https://www.proteinatlas.org/ENSG00000100351-GRAP2: "GRB2 related adaptor protein 2" and "Interacts with SLP-76 to regulate NF-AT activation. Binds to tyrosine-phosphorylated shc"
        'MYL9', # https://www.proteinatlas.org/ENSG00000101335-MYL9: "Myosin light chain 9" and "Myosin regulatory subunit that plays an important role in regulation of both smooth muscle and nonmuscle cell contractile activity via its phosphorylation. Implicated in cytokinesis, receptor capping, and cell locomotion".
        'MYLK', # https://www.proteinatlas.org/ENSG00000065534-MYLK: "Myosin light chain kinase" and "Calcium/calmodulin-dependent myosin light chain kinase implicated in smooth muscle contraction via phosphorylation of myosin light chains (MLC). Also regulates actin-myosin interaction through a non-kinase activity"
        'MAP3K7CL', # https://www.proteinatlas.org/ENSG00000156265-MAP3K7CL: "MAP3K7 C-terminal like" and ""
        'LAT', # https://www.proteinatlas.org/ENSG00000213658-LAT: "Linker for activation of T cells" and "Required for TCR (T-cell antigen receptor)- and pre-TCR-mediated signaling, both in mature T-cells and during their development"
        'MARCH2', # https://www.proteinatlas.org/ENSG00000099785-MARCHF2: "Membrane associated ring-CH-type finger 2" and "E3 ubiquitin-protein ligase that may mediate ubiquitination of TFRC and CD86, and promote their subsequent endocytosis and sorting to lysosomes via multivesicular bodies"
        'ITGA2B', # https://www.proteinatlas.org/ENSG00000005961-ITGA2B: "Integrin subunit alpha 2b" and "Integrin alpha-IIb/beta-3 is a receptor for fibronectin, fibrinogen, plasminogen, prothrombin, thrombospondin and vitronectin. It recognizes the sequence R-G-D in a wide array of ligands. It recognizes the sequence H-H-L-G-G-G-A-K-Q-A-G-D-V in fibrinogen gamma chain. Following activation integrin alpha-IIb/beta-3 brings about platelet/platelet interaction through binding of soluble fibrinogen. This step leads to rapid platelet aggregation which physically plugs ruptured endothelial cell surface".
        'SPARC', # https://www.proteinatlas.org/ENSG00000113140-SPARC: "Secreted protein acidic and cysteine rich" and "Appears to regulate cell growth through interactions with the extracellular matrix and cytokines".
        'PGRMC1', # https://www.proteinatlas.org/ENSG00000101856-PGRMC1: "Progesterone receptor membrane component 1" and "Component of a progesterone-binding protein complex 1. Binds progesterone 2. Has many reported cellular functions (heme homeostasis, interaction with CYPs). Required for the maintenance of uterine histoarchitecture and normal female reproductive lifespan (By similarity). Intracellular heme chaperone. Regulates heme synthesis via interactions with FECH and acts as a heme donor for at least some hemoproteins"
        'PTPN18', # https://www.proteinatlas.org/ENSG00000072135-PTPN18: "Protein tyrosine phosphatase non-receptor type 18" and "Differentially dephosphorylate autophosphorylated tyrosine kinases which are known to be overexpressed in tumor tissues".
        'CCND3', # https://www.proteinatlas.org/ENSG00000112576-CCND3: "Cyclin D3" and "Regulatory component of the cyclin D3-CDK4 (DC) complex that phosphorylates and inhibits members of the retinoblastoma (RB) protein family including RB1 and regulates the cell-cycle during G(1)/S transition".
        'STOM', # https://www.proteinatlas.org/ENSG00000148175-STOM: "Stomatin" and "Regulates ion channel activity and transmembrane ion transport. Regulates ASIC2 and ASIC3 channel activity"
        'RSU1', # https://www.proteinatlas.org/ENSG00000148484-RSU1: "Ras suppressor protein 1" and "Potentially plays a role in the Ras signal transduction pathway. Capable of suppressing v-Ras transformation in vitro".
        'TUBA1C', # https://www.proteinatlas.org/ENSG00000167553-TUBA1C: "Tubulin alpha 1c" and "Tubulin is the major constituent of microtubules. It binds two moles of GTP, one at an exchangeable site on the beta chain and one at a non-exchangeable site on the alpha chain".
        'FYB1', # https://www.proteinatlas.org/ENSG00000082074-FYB1: "FYN binding protein 1" and "Acts as an adapter protein of the FYN and LCP2 signaling cascades in T-cells (By similarity). May play a role in linking T-cell signaling to remodeling of the actin cytoskeleton 1, 2. Modulates the expression of IL2 (By similarity). Involved in platelet activation (By similarity). Prevents the degradation of SKAP1 and SKAP2 3. May be involved in high affinity immunoglobulin epsilon receptor signaling in mast cells".
        'FERMT3', # https://www.proteinatlas.org/ENSG00000149781-FERMT3: "Fermitin family member 3" and "Plays a central role in cell adhesion in hematopoietic cells 1, 2. Acts by activating the integrin beta-1-3 (ITGB1, ITGB2 and ITGB3) (By similarity). Required for integrin-mediated platelet adhesion and leukocyte adhesion to endothelial cells 3. Required for activation of integrin beta-2 (ITGB2) in polymorphonuclear granulocytes (PMNs) (By similarity)".
        'CNST', # https://www.proteinatlas.org/ENSG00000162852-CNST: "Consortin, connexin sorting protein" and "Required for targeting of connexins to the plasma membrane".
        'SNAP23', # https://www.proteinatlas.org/ENSG00000092531-SNAP23: "Synaptosome associated protein 23" and "Essential component of the high affinity receptor for the general membrane fusion machinery and an important regulator of transport vesicle docking and fusion".
        'ETFA', # https://www.proteinatlas.org/ENSG00000140374-ETFA: "Electron transfer flavoprotein subunit alpha" and "Heterodimeric electron transfer flavoprotein that accepts electrons from several mitochondrial dehydrogenases, including acyl-CoA dehydrogenases, glutaryl-CoA and sarcosine dehydrogenase 1, 2, 3, 4, 5. It transfers the electrons to the main mitochondrial respiratory chain via ETF-ubiquinone oxidoreductase (ETF dehydrogenase) 6. Required for normal mitochondrial fatty acid oxidation and normal amino acid metabolism".
        'SVIP', # https://www.proteinatlas.org/ENSG00000198168-SVIP: "Small VCP interacting protein" and "Endoplasmic reticulum-associated degradation (ERAD) is the pathway by which misfolded proteins in the endoplasmic reticulum are targeted to the proteasome for degradation. Multiple specialized proteins interact with one another during ERAD to complete this process. The protein encoded by this gene is an inhibitor of ERAD, functioning to disrupt the interaction of these protein components. This downregulation of ERAD may be needed to protect the cell from overactive protein degradation"
        'PDLIM1', # https://www.proteinatlas.org/ENSG00000107438-PDLIM1: "PDZ and LIM domain 1" and "Cytoskeletal protein that may act as an adapter that brings other proteins (like kinases) to the cytoskeleton 1. Involved in assembly, disassembly and directioning of stress fibers in fibroblasts. Required for the localization of ACTN1 and PALLD to stress fibers. Required for cell migration and in maintaining cell polarity of fibroblasts (By similarity)"
        'BIN2', # https://www.proteinatlas.org/ENSG00000110934-BIN2: "Bridging integrator 2" and "Promotes cell motility and migration, probably via its interaction with the cell membrane and with podosome proteins that mediate interaction with the cytoskeleton. Modulates membrane curvature and mediates membrane tubulation. Plays a role in podosome formation. Inhibits phagocytosis"
        'PIP4K2A', # https://www.proteinatlas.org/ENSG00000150867-PIP4K2A: "Phosphatidylinositol-5-phosphate 4-kinase type 2 alpha" and "Catalyzes the phosphorylation of phosphatidylinositol 5-phosphate (PtdIns5P) on the fourth hydroxyl of the myo-inositol ring, to form phosphatidylinositol 4,5-bisphosphate (PtdIns(4,5)P2) 1, 2. Has both ATP- and GTP-dependent kinase activities".
        'HACD4', # https://www.proteinatlas.org/ENSG00000188921-HACD4: "3-hydroxyacyl-CoA dehydratase 4" and "Catalyzes the third of the four reactions of the long-chain fatty acids elongation cycle. This endoplasmic reticulum-bound enzymatic process, allows the addition of two carbons to the chain of long- and very long-chain fatty acids/VLCFAs per cycle".
        'PKM', # https://www.proteinatlas.org/ENSG00000067225-PKM: "Pyruvate kinase M1/2" and "Catalyzes the final rate-limiting step of glycolysis by mediating the transfer of a phosphoryl group from phosphoenolpyruvate (PEP) to ADP, generating ATP".
        'PDCD10', # https://www.proteinatlas.org/ENSG00000114209-PDCD10: "Programmed cell death 10" and "Promotes cell proliferation. Modulates apoptotic pathways. Increases mitogen-activated protein kinase activity and STK26 activity 1. Important for cell migration, and for normal structure and assembly of the Golgi complex 2. Important for KDR/VEGFR2 signaling. Increases the stability of KDR/VEGFR2 and prevents its breakdown. Required for normal cardiovascular development. Required for normal angiogenesis, vasculogenesis and hematopoiesis during embryonic development (By similarity)".
        'VDAC3', # https://www.proteinatlas.org/ENSG00000078668-VDAC3: "Voltage dependent anion channel 3" and "Forms a channel through the mitochondrial outer membrane that allows diffusion of small hydrophilic molecules".
        'TLK1', # https://www.proteinatlas.org/ENSG00000198586-TLK1: "Tousled like kinase 1" and "Rapidly and transiently inhibited by phosphorylation following the generation of DNA double-stranded breaks during S-phase. This is cell cycle checkpoint and ATM-pathway dependent and appears to regulate processes involved in chromatin assembly".
        'BEX3', # https://www.proteinatlas.org/ENSG00000166681-BEX3: "Brain expressed X-linked 3" and "May be a signaling adapter molecule involved in p75NTR-mediated apoptosis induced by NGF. Plays a role in zinc-triggered neuronal death (By similarity). May play an important role in the pathogenesis of neurogenetic diseases".
        'ITGB1', # https://www.proteinatlas.org/ENSG00000150093-ITGB1: "Integrin subunit beta 1" and "Integrins are heterodimeric proteins made up of alpha and beta subunits. At least 18 alpha and 8 beta subunits have been described in mammals. Integrin family members are membrane receptors involved in cell adhesion and recognition in a variety of processes including embryogenesis, hemostasis, tissue repair, immune response and metastatic diffusion of tumor cells. This gene encodes a beta subunit"
    ],
}
