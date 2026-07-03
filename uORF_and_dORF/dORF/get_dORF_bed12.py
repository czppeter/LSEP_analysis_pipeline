
import re

f3 = open("all_dORF.bed", "w")

utr = {}


with open("utr_3_pos2.txt") as f:
    for line in f:
        tid, blocks, strand = line.strip().split("\t")

        block_list = []

        for b in blocks.split(";"):
            s,e = map(int,b.split("-"))
            block_list.append((s,e))

        utr[tid] = (block_list,strand)

with open("all_dORF_region2.txt") as f:

    for line in f:

        tid, region = line.strip().split("\t")

        orf_start, orf_end = map(int,region.split("-"))
        #考虑到bed12 的格式，在这里加一
        #orf_end = orf_end + 1

        blocks,strand = utr[tid]

        current = 1

        genomic_segments = []

        for bstart,bend in blocks:
            #处理bed文件时已经减去1了
            blen = bend - bstart 

            utr_start = current
            utr_end   = current + blen - 1

            ov_start = max(orf_start,utr_start)
            ov_end   = min(orf_end,utr_end)

            if ov_start <= ov_end:

                offset1 = ov_start - utr_start
                offset2 = ov_end - utr_start

                if strand == "+":
                    #bstart 在前期处理时已近减去了1
                    gstart = bstart + offset1
                    gend   = bstart + offset2 + 1

                else:

                    gstart = bend - offset2 - 1
                    gend   = bend - offset1

                    if gstart > gend:
                        gstart,gend = gend,gstart

                genomic_segments.append((gstart,gend))

            current += blen

        chromStart=min(x[0] for x in genomic_segments)
        chromEnd=max(x[1] for x in genomic_segments)

        blockSizes=[]
        blockStarts=[]

        genomic_segments.sort(key=lambda x: x[0])

        for s,e in genomic_segments:
            blockSizes.append(str(e-s))
            blockStarts.append(str(s-chromStart))
        
        # for i in blockStarts:
        #     if len(blockStarts) == 1:
        #         blockStarts[0] = "0"
        #         break
        #     else:
        #         if i == blockStarts[0]:
        #             blockStarts[0] = "0"
        #         else:
        #             blockStarts[blockStarts.index(i)] = str(int(i) + 1)

        m = re.search(r'LOC_Os(\d+)g', tid)

        chrom = "Chr" + str(int(m.group(1)))


        print(
            chrom,
            chromStart,
            chromEnd,
            tid,
            0,
            strand,
            chromStart,
            chromEnd,
            0,
            len(genomic_segments),
            ",".join(blockSizes),
            ",".join(blockStarts),
            sep="\t"
        )
        bed12 = [chrom,
            str(chromStart),
            str(chromEnd),
            tid,
            "0",
            strand,
            str(chromStart),
            str(chromEnd),
            "0",
            str(len(genomic_segments)),
            ",".join(blockSizes),
            ",".join(blockStarts)]
        

        f3.write("\t".join(bed12) + "\n")
    f3.close


