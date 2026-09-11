# Bundled antibody germline references

Source: [IMGT/GENE-DB](https://www.imgt.org/download/GENE-DB/), release
**202636-7**, downloaded 10 September 2026. IMGT data and metadata are licensed
under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) under the
[IMGT terms retrieved that day](https://www.imgt.org/about/termsofuse.php).
Attribution: IMGT®, the international ImMunoGeneTics information system®,
Montpellier, France. The data retain this license independently of Arpeggia's
source-code license.

The [gapped amino-acid source](https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F%2BORF%2BinframeP)
has SHA-256
`3cb6b0b8cb8940b3b2a9b105771a6a74aa67c06e3ca39eaea0d2030c90e7efd0`.
The subset retains functional IGHV/IGKV/IGLV and IGHJ/IGKJ/IGLJ records for human,
mouse, alpaca, rat and rabbit, including strain/subspecies names. Bracketed/parenthesized
functional labels are included; stop-containing sequences are excluded.
Original headers, partial sequences, ambiguity symbols and IMGT gaps are retained.
There are 1,603 V and 91 J records. Human, mouse, rat and rabbit cover H/K/L;
alpaca references cover heavy chains only. Rat contributes 268 V / 13 J records
and rabbit 123 V / 20 J records.

The subset SHA-256 is
`e87e18cdfae839957b454edb47445d024b95d7eb2ac1b4b655796cad924d8698`.
Regenerate from the pinned download using
[prepare.py](https://github.com/y1zhou/arpeggia/blob/master/data/germlines/prepare.py):

```bash
python3 data/germlines/prepare.py imgt-source.fasta > data/germlines/imgt-202636-7.fasta
```

Matching uses V sequence through IMGT 104 and the J amino-acid sequence. Leading
padding and truncated ends represent unavailable reference coverage. In
particular, alpaca IGHJ5*01 ends before IMGT 128 despite lacking a partial-record
flag; completeness flags alone are insufficient for terminal imputation.
