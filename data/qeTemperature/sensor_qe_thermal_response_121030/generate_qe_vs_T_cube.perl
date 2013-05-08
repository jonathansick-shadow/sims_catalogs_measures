#!/usr/bin/perl -w
# script to repeatedly call ml_design with appropriate parameters, determined
# forensically from vendor provided curves - and output them for use by 
# the Calibration team. This script can be run to regenerate this data.
$cmd="ml_design -m0 vacuum -m1 SiO2 ".
    "-l MgF2 94.47 -l Ta2O5 109.53  -l Si 1.78 -l Si 100000 ".
    "-aa 14.1 23.6 -nw 9001 300 1200 -a";
$wave_field=3;
$deadabs_field=13;
$deplabs_field=15;

for ($temperature=173-5;$temperature<=173+15;$temperature+=0.1) {
    $result=`$cmd -T $temperature`;
    $outfile=sprintf("sensorQE_120227_T%3.1f.txt",$temperature);
    open(G,">$outfile") || die;
    foreach $line ( split("\n",$result) ) {
	@contents=split(' ',$line);
	printf G "%6.2f %8.5f\n",$contents[$wave_field],$contents[$deplabs_field]+0.5*$contents[$deadabs_field];
    }
    close(G);
}
