#DATABASE phreeqd.dat
SOLUTION_MASTER_SPECIES
 A A 0 A 1
SOLUTION_SPECIES
 A = A; log_k 0
 dw 0.9e-9
# H2O + 0.01e- = H2O-0.01; -log_k -12
SOLUTION 0
 pH 7.0 charge
 Cl 1
 Fe 1
 A 1
SOLUTION 1-11
 pH 7 charge
END
PRINT; -reset false
SELECTED_OUTPUT
 -file i3.prn
 -reset false
# -high_prec
USER_PUNCH
 -head dist pH pe O2 A Fe2 Fe3 water
 -start
 1 if step_no < 1 then goto 20
 10 punch dist, -la("H+"), -la("e-"), mol("O2")*1e3, tot("A")*1e3, tot("Fe(2)")*1e3, tot("Fe(3)")*1e3, tot("water") - 1.0
 20 END
 -end
TRANSPORT
 -cells 10
 -shifts 5 0
 -punch_fr 5
 -time 6.6e8
 -bcon 1 1
 -multi_d true 0.9e-9 1.0 0.0 1.0 #   D_w, porosity, porosity_limit, exponent
END
SOLUTION 0
 pH 7.0 charge
 Cl 1
 Fe 1
 A 1
SOLUTION 1-11
 pH 7 charge
TRANSPORT
 -diffc 0.9e-9
 -multi_d false
END