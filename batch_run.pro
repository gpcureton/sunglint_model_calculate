;main,'/data/geoffc/IDL/results/GlintSim_4096_512_4_20000_phillips_3_0_0_nl.hdf','GlintModel_4096_512_4_20000_phillips_3_0_0_nl.hdf',0

;main,'/data/geoffc/IDL/results/GlintSim_4096_512_2_20000_phillips_3_0_0.hdf',   'GlintModel_4096_512_2_20000_phillips_3_0_0.hdf',0
;main,'/data/geoffc/IDL/results/GlintSim_4096_512_2_20000_phillips_3_0_0_nl.hdf','GlintModel_4096_512_2_20000_phillips_3_0_0_nl.hdf',0
;main,'/data/geoffc/IDL/results/GlintSim_4096_512_2_20000_phillips_4_0_0.hdf',   'GlintModel_4096_512_2_20000_phillips_4_0_0.hdf',0
;main,'/data/geoffc/IDL/results/GlintSim_4096_512_2_20000_phillips_4_0_0_nl.hdf','GlintModel_4096_512_2_20000_phillips_4_0_0_nl.hdf',0

;main,'/data/geoffc/IDL/results/GlintSim_4096_512_4_20000_phillips_3_0_0.hdf',   'GlintModel_4096_512_4_20000_phillips_3_0_0.hdf',0
;main,'/data/geoffc/IDL/results/GlintSim_4096_512_4_20000_phillips_3_0_0_nl.hdf','GlintModel_4096_512_4_20000_phillips_3_0_0_nl.hdf',0
;main,'/data/geoffc/IDL/results/GlintSim_4096_512_4_20000_phillips_4_0_0.hdf',   'GlintModel_4096_512_4_20000_phillips_4_0_0.hdf',0
;main,'/data/geoffc/IDL/results/GlintSim_4096_512_4_20000_phillips_4_0_0_nl.hdf','GlintModel_4096_512_4_20000_phillips_4_0_0_nl.hdf',0

;main,'/data/geoffc/IDL/results/GlintSim_4096_512_2_160000_phillips_3_0_0.hdf',   'GlintModel_4096_512_2_160000_phillips_3_0_0.hdf',0
;main,'/data/geoffc/IDL/results/GlintSim_4096_512_2_160000_phillips_3_0_0_nl.hdf','GlintModel_4096_512_2_160000_phillips_3_0_0_nl.hdf',0

;main,'/data/geoffc/IDL/results/GlintSim_4096_512_2_160000_phillips_4_0_0.hdf',     '/data/geoffc/IDL/results/GlintModel_4096_512_2_160000_phillips_4_0_0.hdf',0
main,'/data/geoffc/IDL/results/GlintSim_4096_512_2_160000_phillips_4_0_0_nl.hdf',     '/data/geoffc/IDL/results/GlintModel_4096_512_2_160000_phillips_4_0_0_nl.hdf',0
;main,'/data/geoffc/IDL/results/GlintSim_4096_512_4_160000_phillips_3_0_0.hdf',     '/data/geoffc/IDL/results/GlintModel_4096_512_4_160000_phillips_3_0_0.hdf',0
;main,'/data/geoffc/IDL/results/GlintSim_4096_512_4_160000_phillips_3_0_0_nl.hdf',  '/data/geoffc/IDL/results/GlintModel_4096_512_4_160000_phillips_3_0_0_nl.hdf',0


;inFile = 'GlintSim_4096_512_2_20000_phillips_3_0_0.hdf'
;inFile = 'GlintSim_4096_512_4_20000_phillips_3_0_0.hdf'

;;inFile = 'GlintSim_4096_512_2_20000_phillips_3_0_0_nl.hdf'
;inFile = 'GlintSim_4096_512_4_20000_phillips_3_0_0_nl.hdf'

;inPath = '../results/exponent_three/bareB2_L.bareB3_L.M2_L-sub/'
;outPath =               './output/exponent_three/bareB2_L.bareB3_L.M2_L-sub/'

;inPath = '/data/geoffc/IDL/results/'
;outPath = '/data/geoffc/IDL/results/'

;inFullPath = STRCOMPRESS(inPath+inFile)
;outFile= STRCOMPRESS('GlintModel'+strsplit(inFile,'GlintSim',/REGEX,/EXTRACT))
;outFullPath = STRCOMPRESS(outPath+outFile)
;outFullPath = STRING(STRCOMPRESS(outFile))

;main,inFullPath,outFullPath,0
;main,inFullPath,outFullPath,5
;main,inFullPath,outFullPath,13
