%path = @runpath
cd %path
load Return.wf1

pageselect data

smpl @all
equation garch_vol.arch(tdist) stock c ar(1)

garch_vol.makegarch stock_vol

freeze(volatiltiy) stock_vol
volatiltiy.save(t=csv) stock_vol.csv

