# HLTrigger Study code

Info from the trigger tutorial talk: https://indico.cern.ch/event/375107/contribution/6/material/slides/0.pdf

## Setup
```bash
cmsrel CMSSW_7_4_0_pre9
cd CMSSW_7_4_0_pre9/src
cmsenv
git clone https://github.com/DESY-CMS-SUS/TriggerStudy.git
scram b -j4
```

## Content Description
You’ll	find	(under	TriggerStudy/SimpleHLTAnalyzer)

* in	plugins:	a	simple	EDAnalyzer	to	analyze	the	trigger	results	and	make	a	TTree	(more	in	backup	on	how	this	is	done)	
* in	test:	a	config	which	will	run	the	HLT	process	including	our	trigger	and	auxiliary	trigger	on	RAW	input,	the	EDAnalyzer	is	added	to	the	process	to	produce	the	TTree	(more	in	backup	on	how	to	prepare	this)	
* in	macros:	a	ROOT	macro	to	analyze	the	TTree,	compute	the	efficiencies,	
and	make	plots

## Running the code:
```bash
cd TriggerStudy/SimpleHLTAnalyzer/macros
root -l -q -b makeTrigEffPlots.C+
```
