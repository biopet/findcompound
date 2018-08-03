organization := "com.github.biopet"
organizationName := "Biopet"

startYear := Some(2018)

name := "FindCompound"
biopetUrlName := "findcompound"

biopetIsTool := true

mainClass in assembly := Some(
  s"nl.biopet.tools.${name.value.toLowerCase()}.${name.value}")

developers += Developer(id = "ffinfo",
                        name = "Peter van 't Hof",
                        email = "pjrvanthof@gmail.com",
                        url = url("https://github.com/ffinfo"))

scalaVersion := "2.11.12"

libraryDependencies += "com.github.biopet" %% "tool-utils" % "0.4"
libraryDependencies += "com.github.biopet" %% "ngs-utils" % "0.4"
libraryDependencies += "com.github.biopet" %% "tool-test-utils" % "0.2.2" % Test
libraryDependencies += "com.github.broadinstitute" % "picard" % "2.18.7"
