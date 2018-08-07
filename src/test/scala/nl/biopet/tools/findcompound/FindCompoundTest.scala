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

import java.io.File

import nl.biopet.utils.io.getLinesFromFile
import nl.biopet.utils.test.tools.ToolTest
import org.testng.annotations.Test

class FindCompoundTest extends ToolTest[Args] {
  def toolCommand: FindCompound.type = FindCompound
  @Test
  def testNoArgs(): Unit = {
    intercept[IllegalArgumentException] {
      FindCompound.main(Array())
    }
  }

  @Test
  def testDefault(): Unit = {
    val outputDir = File.createTempFile("test.", ".out")
    outputDir.delete()
    outputDir.mkdir()
    outputDir.deleteOnExit()
    FindCompound.main(
      Array("-i",
            resourcePath("/multi.vcf.gz"),
            "-o",
            outputDir.getAbsolutePath,
            "-R",
            resourcePath("/fake_chrQ.fa"),
            "-r",
            resourcePath("/chrQ.refflat")))

    new File(outputDir, "exon.counts") should exist
    new File(outputDir, "intron.counts") should exist
    new File(outputDir, "total.counts") should exist

    val header = "#Gene\tCompound\tHomVar\t" +
      "Sample_1-het\tSample_1-homVar\tSample_1-homRef\t" +
      "Sample_2-het\tSample_2-homVar\tSample_2-homRef\t" +
      "Sample_3-het\tSample_3-homVar\tSample_3-homRef"

    getLinesFromFile(new File(outputDir, "exon.counts"))
      .mkString("\n") shouldBe List(
      header,
      "geneA\t1\t0\t1\t1\t0\t1\t0\t1\t0\t0\t2",
      "geneB\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0").mkString("\n")
    getLinesFromFile(new File(outputDir, "intron.counts"))
      .mkString("\n") shouldBe List(
      header,
      "geneA\t1\t0\t1\t0\t0\t0\t0\t1\t0\t1\t0",
      "geneB\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0").mkString("\n")
    getLinesFromFile(new File(outputDir, "total.counts"))
      .mkString("\n") shouldBe List(
      header,
      "geneA\t2\t0\t2\t1\t0\t1\t0\t2\t0\t1\t2",
      "geneB\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0").mkString("\n")
  }

  @Test
  def testPedigree(): Unit = {
    val outputDir = File.createTempFile("test.", ".out")
    outputDir.delete()
    outputDir.mkdir()
    outputDir.deleteOnExit()
    FindCompound.main(
      Array(
        "-i",
        resourcePath("/multi.vcf.gz"),
        "-o",
        outputDir.getAbsolutePath,
        "-R",
        resourcePath("/fake_chrQ.fa"),
        "-r",
        resourcePath("/chrQ.refflat"),
        "-p",
        resourcePath("/pedFile.ped")
      ))

    new File(outputDir, "exon.counts") should exist
    new File(outputDir, "intron.counts") should exist
    new File(outputDir, "total.counts") should exist
  }

  @Test
  def testWrongPedigree(): Unit = {
    val outputDir = File.createTempFile("test.", ".out")
    outputDir.delete()
    outputDir.mkdir()
    outputDir.deleteOnExit()
    intercept[IllegalArgumentException] {
      FindCompound.main(
        Array(
          "-i",
          resourcePath("/multi.vcf.gz"),
          "-o",
          outputDir.getAbsolutePath,
          "-R",
          resourcePath("/fake_chrQ.fa"),
          "-r",
          resourcePath("/chrQ.refflat"),
          "-p",
          resourcePath("/wrong.ped")
        ))
    }.getMessage shouldBe "Samples in vcf file not found in ped file: Sample_1, Sample_2, Sample_3"
  }
}
