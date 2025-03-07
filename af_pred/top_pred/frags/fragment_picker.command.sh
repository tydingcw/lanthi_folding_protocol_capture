#!/bin/bash

PEP=$1

/dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin/fragment_picker.linuxgccrelease \
	@ fragment_picker_quota.options \
	-in:file:checkpoint ${PEP}.checkpoint \
	-in:file:fasta ${PEP}.fasta \
	-frags::ss_pred ${PEP}.af_ss jufo \
	-out:file:frag_prefix ${PEP}_frags \
	-frags::describe_fragments ${PEP}_frags.fsc
