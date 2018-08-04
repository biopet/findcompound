package nl.biopet.tools

import picard.annotation.Gene

package object findcompound {
  case class VariantCounts(var het: Int = 0,
                           var homVar: Int = 0,
                           var homRef: Int = 0) {
    def counts: List[Int] = List(het, homVar, homRef)
  }
  case class VariantTypes(exon: VariantCounts = VariantCounts(),
                          intron: VariantCounts = VariantCounts()) {
    def exonCounts: List[Int] = exon.counts
    def intronCounts: List[Int] = intron.counts
    def totalCounts: List[Int] = exon.counts.zip(intron.counts).map{ case (e, i) => e + i}
  }

  case class GeneResults(gene: Gene, samples: IndexedSeq[VariantTypes]) {
    def exonHomVarCount: Int = samples.count(_.exon.homVar > 0)
    def exonCompoundCount: Int = samples.count(x => x.exon.homVar == 0 && x.intron.het >= 2)

    def intronHomVarCount: Int = samples.count(_.intron.homVar > 0)
    def intronCompoundCount: Int = samples.count(x => x.intron.homVar == 0 && x.intron.het >= 2)

    def totalHomVarCount: Int = samples.count(x => x.exon.homVar + x.intron.homVar > 0)
    def totalCompoundCount: Int = samples.count(x => x.exon.homVar + x.intron.homVar == 0 && x.exon.het + x.intron.het >= 2)
  }
}
