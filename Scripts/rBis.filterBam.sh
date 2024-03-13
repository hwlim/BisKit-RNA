#!/usr/bin/env bash

## exit if any error occurs
set -e

alignDir=$1
bam=$2

## get relevant reads for miRNA m5C candidate calling
echo "Get SAM header..."
samtools view -H $bam > $alignDir/header.sam

### split to forward strand
samtools view -b -F 16 $bam > $alignDir/align.plus.bam
samtools index $alignDir/align.plus.bam

### split to minus strand
samtools view -b -f 16 $bam > $alignDir/align.minus.bam
samtools index $alignDir/align.minus.bam


echo "Getting uniquely aligned forward reads..."
samtools view $alignDir/align.plus.bam | grep -w "NH:i:1" > $alignDir/align.plus.uniq.sam || true
cat $alignDir/header.sam $alignDir/align.plus.uniq.sam | samtools view -b - > $alignDir/align.plus.uniq.bam
samtools index $alignDir/align.plus.uniq.bam

echo "Getting uniquely aligned reverse reads..."
samtools view $alignDir/align.minus.bam | grep -w "NH:i:1" > $alignDir/align.minus.uniq.sam || true
cat $alignDir/header.sam $alignDir/align.minus.uniq.sam | samtools view -b - > $alignDir/align.minus.uniq.bam

echo "Getting multi-aligned forward reads..."
samtools view $alignDir/align.plus.bam | grep -wv "NH:i:1" > $alignDir/align.plus.multi.sam || true
cat $alignDir/header.sam $alignDir/align.plus.multi.sam | samtools view -b - > $alignDir/align.plus.multi.bam

echo "Getting multi-aligned reverse reads..."
samtools view $alignDir/align.minus.bam | grep -wv "NH:i:1" > $alignDir/align.minus.multi.sam || true
cat $alignDir/header.sam $alignDir/align.minus.multi.sam | samtools view -b - > $alignDir/align.minus.multi.bam

## Among forward/plus multi aligned reads, get reads aligned exactly once to the forward strand
samtools view $alignDir/align.plus.multi.bam | cut -f1 | sort | uniq -u > $alignDir/forwardUniqReadNames.txt

## use resulting txt file to get reads aligned to the forward strand exactly once
samtools view -N $alignDir/forwardUniqReadNames.txt -o $alignDir/align.plus.exactly_once.bam $alignDir/align.plus.multi.bam

## Among reverse/minus multi aligned reads, get reads aligned exactly once to the reverse strand
samtools view $alignDir/align.minus.multi.bam | cut -f1 | sort | uniq -u > $alignDir/reverseUniqReadNames.txt

## use resulting txt file to get reads aligned to the reverse strand exactly once
samtools view -N $alignDir/reverseUniqReadNames.txt -o $alignDir/align.minus.exactly_once.bam $alignDir/align.minus.multi.bam

## merge resulting bam and uniq.bam
samtools merge -f $alignDir/merged.bam $alignDir/align.plus.uniq.bam $alignDir/align.plus.exactly_once.bam
samtools sort $alignDir/merged.bam > $alignDir/align.plus.filtered.bam
samtools index $alignDir/align.plus.filtered.bam

rm $alignDir/align.plus.uniq.sam
rm $alignDir/align.plus.uniq.bam
rm $alignDir/align.plus.uniq.bam.bai
rm $alignDir/align.minus.uniq.sam
rm $alignDir/align.minus.uniq.bam
rm $alignDir/align.plus.exactly_once.bam
rm $alignDir/align.minus.exactly_once.bam
rm $alignDir/align.plus.multi.sam
rm $alignDir/align.plus.multi.bam
rm $alignDir/align.minus.multi.sam
rm $alignDir/align.minus.multi.bam
rm $alignDir/forwardUniqReadNames.txt
rm $alignDir/reverseUniqReadNames.txt
rm $alignDir/merged.bam
rm $alignDir/header.sam
