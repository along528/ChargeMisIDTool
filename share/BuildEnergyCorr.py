from ROOT import TH1F
from array import array

#these are taken from figure 13 in the same-sign WW support note (ATL-COM-PHYS-2013-990)
#https://cds.cern.ch/record/1561731/files/ATL-COM-PHYS-2013-990.pdf


etaBins = [0.,0.5,1.,1.38,1.9,2.2,2.7]
hEnergyCorrection = TH1F("hEnergyCorrection","hEnergyCorrection; |#eta|;Energy Bias",len(etaBins)-1,array('d',etaBins))
hEnergySmearing = TH1F("hEnergySmearing","hEnergySmearing; |#eta|;Energy Smearing",len(etaBins)-1,array('d',etaBins))


hEnergyCorrection.SetBinContent(1,1.6)
hEnergyCorrection.SetBinContent(2,2.2)
hEnergyCorrection.SetBinContent(3,3.05)
hEnergyCorrection.SetBinContent(4,4.5)
hEnergyCorrection.SetBinContent(5,4.9)
hEnergyCorrection.SetBinContent(6,3.3)

hEnergySmearing.SetBinContent(1,1.75)
hEnergySmearing.SetBinContent(2,2.1)
hEnergySmearing.SetBinContent(3,2.6)
hEnergySmearing.SetBinContent(4,3.45)
hEnergySmearing.SetBinContent(5,3.55)
hEnergySmearing.SetBinContent(6,3.75)


from RootOutput import RootOutput

output = RootOutput("EnergyCorrections.root")
output.write(hEnergyCorrection)
output.write(hEnergySmearing)
