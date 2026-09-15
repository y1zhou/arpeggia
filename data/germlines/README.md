# Bundled antibody germline references

IMGT data and metadata are licensed
under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) under the
[IMGT terms](https://www.imgt.org/about/termsofuse.php).
Attribution: IMGT®, the international ImMunoGeneTics information system®,
Montpellier, France. The data retain this license independently of Arpeggia's
source-code license.

## IMGT/GENE-DB bulk subset

Source: [IMGT/GENE-DB](https://www.imgt.org/download/GENE-DB/), release
**202636-7**, downloaded 10 September 2026.

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

## Llama protein-display supplement

The bulk snapshot and GENE-DB exports queried on 15 September 2026 contained
no `Lama glama` records. The separate
[llama FASTA](https://github.com/y1zhou/arpeggia/blob/master/data/germlines/imgt-llama-20260915.fasta)
comes from IMGT's
[IGHV](https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?group=IGHV&latin=Lama+glama&species=llama)
and [IGHJ](https://www.imgt.org/IMGTrepertoire/Proteins/proteinDisplays.php?group=IGHJ&latin=Lama+glama&species=llama)
protein displays, retrieved 15 September 2026. It retains six functional V
references (AF305944–AF305949) and five J references (AF305952); no light chains
are included. This small historical set is not a comprehensive llama repertoire.
IGHJ5/6 end at IMGT 127, leaving position 128 unavailable for imputation.

[prepare_llama.py](https://github.com/y1zhou/arpeggia/blob/master/data/germlines/prepare_llama.py)
removes HTML and layout spaces, preserves IMGT dots and retains functional rows.
The six header fields are accession, gene/allele, species, functionality, domain
label and `protein-display`. The last field identifies the source type because
these tables do not supply nucleotide spans. It is part of the stable reference ID.
No accession coordinates are inferred. SHA-256 checksums:

| File | SHA-256 |
| --- | --- |
| V display | `e36c191249e7b67ea31da66e848e07ba09605a388c617842514f378d7cd005b7` |
| J display | `ec579819eef533d4031c857de310582d1663982b7d1ab23f9598273edd035e49` |
| Generated FASTA | `016d83d153cbb9a11db498b694932089b440046190c51cec92d527d7f2fa6b10` |

```bash
curl --fail --location 'https://www.imgt.org/3Dstructure-DB/cgi/DomainDisplay-include.cgi?species=Lama+glama&groups=IGHV' -o llama-v.html
curl --fail --location 'https://www.imgt.org/3Dstructure-DB/cgi/DomainDisplay-include.cgi?species=Lama+glama&groups=IGHJ' -o llama-j.html
python3 data/germlines/prepare_llama.py llama-v.html llama-j.html > data/germlines/imgt-llama-20260915.fasta
```

Together, the two files contain 1,609 V and 96 J reference records.
