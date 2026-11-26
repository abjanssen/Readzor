#!/bin/bash
exec micromamba run -p ./envs/mm_fastp fastp "$@"
