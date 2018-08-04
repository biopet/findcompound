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
    def totalCounts: List[Int] =
      exon.counts.zip(intron.counts).map { case (e, i) => e + i }
  }

  case class GeneResults(gene: Gene, samples: IndexedSeq[VariantTypes]) {
    def exonHomVarCount: Int = samples.count(_.exon.homVar > 0)
    def exonCompoundCount: Int =
      samples.count(x => x.exon.homVar == 0 && x.intron.het >= 2)

    def intronHomVarCount: Int = samples.count(_.intron.homVar > 0)
    def intronCompoundCount: Int =
      samples.count(x => x.intron.homVar == 0 && x.intron.het >= 2)

    def totalHomVarCount: Int =
      samples.count(x => x.exon.homVar + x.intron.homVar > 0)
    def totalCompoundCount: Int =
      samples.count(x =>
        x.exon.homVar + x.intron.homVar == 0 && x.exon.het + x.intron.het >= 2)
  }
}
