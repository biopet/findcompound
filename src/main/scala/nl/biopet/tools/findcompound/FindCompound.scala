/*
 * Copyright (c) 2018 Biopet
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package nl.biopet.tools.findcompound

import java.io.{File, PrintWriter}

import htsjdk.variant.vcf.VCFFileReader
import nl.biopet.utils.ngs.intervals.BedRecord
import nl.biopet.utils.ngs.{fasta, vcf}
import nl.biopet.utils.tool.ToolCommand
import picard.annotation.{Gene, GeneAnnotationReader}

import scala.collection.JavaConversions._

object FindCompound extends ToolCommand[Args] {
  def emptyArgs = Args()
  def argsParser = new ArgsParser(this)

  def main(args: Array[String]): Unit = {
    val cmdArgs = cmdArrayToArgs(args)

    logger.info("Start")

    //Sets picard logging level
    htsjdk.samtools.util.Log
      .setGlobalLogLevel(
        htsjdk.samtools.util.Log.LogLevel.valueOf(logger.getLevel.toString))

    logger.info("Reading refflat file")

    val geneReader = GeneAnnotationReader.loadRefFlat(
      cmdArgs.refflatFile,
      fasta.getCachedDict(cmdArgs.referenceFasta))

    val vcfReader = new VCFFileReader(cmdArgs.inputVcfFile, true)

    val samples = vcf.getSampleIds(cmdArgs.inputVcfFile).toIndexedSeq
    val sampleMap = samples.zipWithIndex.toMap

    val results = geneReader.getAll
      .map(createGeneResults(_, sampleMap, vcfReader))
      .toList
      .sortBy(_.gene.getName)

    writeOutput(results,
                samples,
                new File(cmdArgs.outputDir, "exon.counts"),
                _.exonHomVarCount,
                _.exonCompoundCount,
                _.exonCounts)
    writeOutput(results,
                samples,
                new File(cmdArgs.outputDir, "intron.counts"),
                _.intronHomVarCount,
                _.intronCompoundCount,
                _.intronCounts)
    writeOutput(results,
                samples,
                new File(cmdArgs.outputDir, "total.counts"),
                _.totalHomVarCount,
                _.totalCompoundCount,
                _.totalCounts)

    logger.info("Done")
  }

  /** This method reads the vcf file per given gene and create a [[GeneResults]] */
  def createGeneResults(gene: Gene,
                        sampleMap: Map[String, Int],
                        vcfReader: VCFFileReader): GeneResults = {
    val r =
      GeneResults(gene, IndexedSeq.fill(sampleMap.size)(VariantTypes()))
    val variants =
      vcf.loadRegion(vcfReader,
                     BedRecord(gene.getContig, gene.getStart, gene.getEnd))
    val exons = gene.flatMap(_.exons).toList

    variants.foreach { variant =>
      variant.getGenotypes.filter(_.isCalled).foreach { g =>
        val idx = sampleMap(g.getSampleName)
        val exonic = exons.exists(e =>
          e.start <= variant.getStart && e.end >= variant.getEnd)
        g match {
          case _ if g.isHet && exonic    => r.samples(idx).exon.het += 1
          case _ if g.isHet              => r.samples(idx).intron.het += 1
          case _ if g.isHomRef && exonic => r.samples(idx).exon.homRef += 1
          case _ if g.isHomRef           => r.samples(idx).intron.homRef += 1
          case _ if g.isHomVar && exonic => r.samples(idx).exon.homVar += 1
          case _ if g.isHomVar           => r.samples(idx).intron.homVar += 1
          case _                         => throw new IllegalStateException("Please fix this")
        }
      }
    }
    r
  }

  /** This method will write results to a file */
  def writeOutput(results: List[GeneResults],
                  samples: IndexedSeq[String],
                  outputFile: File,
                  homVarCount: GeneResults => Int,
                  compoundCount: GeneResults => Int,
                  sampleCounts: VariantTypes => List[Int]): Unit = {
    val headerLine = (List("#Gene", "Compound", "HomRef") ++ samples.flatMap(
      s => List(s"$s-het", s"$s-homVar", s"$s-homRef"))).mkString("\t")

    val writer = new PrintWriter(outputFile)
    writer.println(headerLine)

    results.foreach { gene =>
      val geneName = gene.gene.getName

      writer.println(
        (List(geneName, homVarCount(gene), compoundCount(gene)) ++ gene.samples
          .flatMap(sampleCounts))
          .mkString("\t"))
    }
    writer.close()
  }

  def descriptionText: String =
    """
      |this tool will count the number of samples that has a homozygous event or a compound homozygous event per gene.
      |It's advised to first filter on annotation scores such as sift and polyphen before using this tool.
    """.stripMargin

  def manualText: String =
    """
      |To use this the reference file should have a dict file next to it.
      |The vcf file should have genotypes and should be in the same contigs names as the reference.
    """.stripMargin

  def exampleText: String =
    s"""
      |Default run with a vcf file:
      |${example("-i",
                 "<input vcf file>",
                 "-o",
                 "<output dir>",
                 "-R",
                 "<reference fasta>",
                 "-r",
                 "<refflat file>")}
      |
    """.stripMargin
}
