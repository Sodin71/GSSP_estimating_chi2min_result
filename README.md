# GSSP_estimating_chi2min_result
This program is for estimating results from the grid-search result file of GSSP

created in 2023-12-09

## Usage
### input file: Chi2_table.dat (GSSP result file)
### output file: ".txt" and ".png" files for each parameter
- The files are only created for "adjusted" parameters in GSSP.

## optional
- The results are also listed in the console/terminal
- The default dimension value is 3. You can change it if you want (check the comment after "if __name__ == "__main__":" line)
- "input file" can change (check the comment after "if __name__ == "__main__":" line)
- This program works for 5 parameters: 'MH', 'Teff', 'log g', 'micro_turbulence', and 'v sin i'. You can select a parameter to work (check the comment after "if __name__ == "__main__":" line).
