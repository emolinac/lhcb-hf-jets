## Description of important scripts
**Note:** All scripts that are run to produce the final results are shown for informational purposes, but in practice, convenient scripts that streamline the analysis are used. See the **analysis workflow** section for details
- **`MCMakeVarTree.C`**: Used to create slimmed down TTrees from the DaVinci output for truth MC. Each entry of the resulting tree corresponds to a B-jet candidate
  ```bash
  # Create truth MC TTrees
  root -l -b -q "MCMakeVarTree.C(-1, 61590)"
  root -l -b -q "MCMakeVarTree.C(-1, 61591)"
  root -l -b -q "MCMakeVarTree.C(-1, 71590)"
  root -l -b -q "MCMakeVarTree.C(-1, 71591)"
  root -l -b -q "MCMakeVarTree.C(-1, 81590)"
  root -l -b -q "MCMakeVarTree.C(-1, 81591)"
  ```
- **`MCSimpleEfficiency.C`**: Processes output of `MCMakeVarTree.C` to obtain the denominator for the efficiency calculation (efficiency = # of matched B-jets / # of truth B-jets)
  ```bash
  root -l -b -q "MCSimpleEff.C(-1, 91599)"
  ```
- **`MCSimpleObservalbes.C`**: Processes output of `MCMakeVarTree.C` to obtain truth level distributions used for the (1) closure tests and (2) comparison to final results
  ```bash
  root -l -b -q "MCMakeObservable.C(-1, 91599)"
  ```
- **`MakeVarTree.C`**: Used to create slimmed down TTrees from the DaVinci output for data and reco MC. Each entry of the resulting tree corresponds to a B-jet candidate
  ```bash
  # Create data TTrees
  root -l -b -q "MakeVarTree.C(-1, 61590, 1)"
  root -l -b -q "MakeVarTree.C(-1, 61591, 1)"
  root -l -b -q "MakeVarTree.C(-1, 71590, 1)"
  root -l -b -q "MakeVarTree.C(-1, 71591, 1)"
  root -l -b -q "MakeVarTree.C(-1, 81590, 1)"
  root -l -b -q "MakeVarTree.C(-1, 81591, 1)"
  
  #Create reco MC TTrees
  root -l -b -q "MakeVarTree.C(-1, 61590, 0)"
  root -l -b -q "MakeVarTree.C(-1, 61591, 0)"
  root -l -b -q "MakeVarTree.C(-1, 71590, 0)"
  root -l -b -q "MakeVarTree.C(-1, 71591, 0)"
  root -l -b -q "MakeVarTree.C(-1, 81590, 0)"
  root -l -b -q "MakeVarTree.C(-1, 81591, 0)"
  ```
- **`MassFit.C`**: Processes output of `MakeVarTree.C` to fit the $B^{\pm}$ invariant mass. Requires running workflow in the `../BjetMisID` directory (see [../BjetMisID](https://gitlab.cern.ch/LHCb-QEE/hfjetfragmentation/-/tree/main/BjetMisID))
  ```bash
  root -l -b -q "MassFit.C(-1, 91599, 0)" # reco MC
  root -l -b -q "MassFit.C(-1, 91599, 1)" # data
  ```
- **`MakeCorrections.C`**: Processes output of `MakeVarTree.C` and `MCSimpleEfficiency.C` for reco MC to obtain (1) purities (purity = # of matched B-jets / # of reco B-jets), (2) response matrices for the unfolding (from matched B-jets), and (3) efficincies (numerator of "# of matched B-jets" combined with denominator from `MCSimpleEfficiencies.C` output)
    - Reads in output of `MassFit.C` to make the relevant mass cut
  ```bash
  root -l -b -q "MakeCorrections.C(-1, 91599)"
  ```
- **`MakeObservable.C`**: Processes output of `MakeVarTree.C` to obtain observable distributions for data and reco MC
    - Reads in histograms for trigger / track / pid corrections
    - Reads in output of `MassFit.C` to make the relevant mass cut and perform the sideband subtraction
  ```bash
  root -l -b -q "MakeObservable.C(-1, 91599, 0)"
  root -l -b -q "MakeObservable.C(-1, 91599, 1)"
  ```
- **`ClosureTest.C`**: Processes output of `MakeObservable.C` for reco MC, `MakeCorrections.C`, and `MCMakeObservable.C` to perform both statistical and shape closure test depending on user arguments 
  ```bash
  # statistical closure test
  root -l -b -q "ClosureTest.C(-1, 91599, 91599, 4, 0)"
  # shape closure test
  root -l -b -q "ClosureTest.C(-1, 91599, 91599, 4, 1)"
  ```
- **`MakeRMReweightFactors.C`**: Creates ratios of $(z, j_{T}, p_{T,\text{jet}})$ distributions in data / reco MC. Used to reweight the response matrix for the prior and shape nonclosure systematic unceratinties
  ```bash
  root -l -b -q MakeRMReweightFactors.C
  ```
- **`PlotFinal.C`**: Processes output of `MakeObservable.C` and `MakeCorrections.C` to apply the correction chain (purity -> unfolding -> efficiency) to observable distributions in data.
  ```bash
  root -l -b -q PlotFinal.C
  ```
- **`Systematics/GetSystematic.C`**: Processes output from `MakeObservable.C` and `MakeCorrections.C` ran with approrpiate user arguments for relevant systematic variations to compute systematic uncertainties for the various sources. 
    - User can provide arguments to calculate uncertainties for the different systematic variations. Preparing these inputs is discussed in a subsection below.
  ```bash
  cd Systematics/
  #JESJER
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)"
  #JetID
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)"
  #TrackEff Up
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)"
  #TrackEff Down
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0)"
  #PIDEff Up
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)"
  #PIDEff Down
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0)"
  #TrigEff Up
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)"
  #TrigEff Down
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0)"
  #RecSelEff
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)"
  #Fit Model
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)"
  #SB Sub Near
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0)"
  #SB Sub Far
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0)"
  #Iter Sys Up
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)"
  #Iter Sys Down
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0)"
  #Prior
  root -l -b -q "GetSystematic.C(-1, 91599, 91599, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)"
  cd ../
  ```
- **`Systematics/GetTotalSys.C`**: Processes output from `Systematics/GetSystematic.C` from **all** relevant systematic variations to obtain the total systematic uncertainties for each observable / pair of observables
  ```bash
  cd Systematics/
  root -l -b -q GetTotalSys.C
  cd ../
  ```
- **`PublicationFigures.C`**: Processes output from `PlotFinal.C`, `Systematics/GetTotalSys.C`, and `MCMakeObservable.C` to obtain final observable distributions with systematic uncertainties and comparison to truth MC
  ```bash
  root -l -b -q PublicationFigures.C
  ```

## Analysis workflow

**Additional System Requirement:** `tmux`
- Terminal multiplexer -- easy solution for running scripts in parallel when possible

In practice, bash scripts are used to help automate / streamline the workflow, and files are processed in parallel when possible to speed up the analysis chain. The following scripts run the full analysis chain described above:

**Note:** This assumes that you have already processed the inputs to the MisID mass fit! 

```bash
# Create truth MC TTrees for all years and magnet polarities
source runMCMakeVarTree.sh
tmux a #(attach to tmux session)
#wait until all scripts finish running...
ctrl+D #(detatch from tmux session)

# Create reco MC TTrees for all years and magnet polarities
source runMakeVarTreeReco.sh
tmux a #(attach to tmux session)
#wait until all scripts finish running...
ctrl+D #(detatch from tmux session)

# Create data TTrees for all years and magnet polarities
source runMakeVarTreeData.sh
tmux a #(attach to tmux session)
#wait until all scripts finish running...
ctrl+D #(detatch from tmux session)

# Process TTrees to obtain observable distributions and corresponding correction histograms
source processTrees.sh

# Statitsical closure test
root -l -b -q "ClosureTest.C(-1, 91599, 91599, 4, 0)"

# Shape closure test
root -l -b -q "ClosureTest.C(-1, 91599, 91599, 4, 1)"

cd Systematics/

# Create TTrees with systematic variations (2016)
source runAllMakeVarTree_2016_Systematics.sh
tmux a #(attach to tmux session)
#wait until all scripts finish running...
ctrl+D #(detatch from tmux session)

# Create TTrees with systematic variations (2017)
source runAllMakeVarTree_2017_Systematics.sh
tmux a #(attach to tmux session)
#wait until all scripts finish running...
ctrl+D #(detatch from tmux session)

# Create TTrees with systematic variations (2018)
source runAllMakeVarTree_2018_Systematics.sh
tmux a #(attach to tmux session)
#wait until all scripts finish running...
ctrl+D #(detatch from tmux session)

# Systematic variations for mass fits
source runAllMassFit_Systematics.sh

# Systematic variations for purities, response, efficiencies
source runAllMakeCorrections_Systematics.sh

# Systematic variations for observable distributions
source runAllMakeObservable_Systematics.sh

# Calcualte relative systematic uncertainty from the systematic variations (1)
source CalculateSystematics1.sh
tmux a #(attach to tmux session)
#wait until all scripts finish running...
ctrl+D #(detatch from tmux session)

# Calcualte relative systematic uncertainty from the systematic variations (2)
source CalculateSystematics2.sh
tmux a #(attach to tmux session)
#wait until all scripts finish running...
ctrl+D #(detatch from tmux session)

# Combine systematic uncertainties to obtain the total relative systematic uncertainty
root -l -b -q GetTotalSys.C

cd ../

# Create final distriubtions
root -l -b -q PublicationFigures.C
```