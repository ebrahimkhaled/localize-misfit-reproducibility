# worked_example.R -- the worked example of the paper's software section, run with the package functions: the SUPPORT model M0
# developed on acute diagnoses and validated on chronic ones. The seeds are those of apps_update.R and size_apps.R, so the
# output must equal the paper's verdict and sizes. Output worked_example.txt.
suppressMessages({ library(splines); devtools::load_all("C:/Users/ebrah/.cursor-tutor/projects/ebgof-release-2.10.0", quiet = TRUE) })
HERE <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/learned-def/option_a"
s <- read.csv("C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/paper_deepgof/medical_example/data/support2.csv",
              na.strings = c("NaN", "NA", ""))
d <- na.omit(data.frame(y = as.integer(s$hospdead), age = s$age, meanbp = s$meanbp, hrt = s$hrt, resp = s$resp, temp = s$temp,
                        crea = s$crea, sod = s$sod, wblc = s$wblc, num.co = s$num.co, male = as.numeric(s$sex == "male"),
                        acute = as.integer(s$dzclass %in% c("ARF/MOSF", "Coma"))))
acute <- d[d$acute == 1, ]; chronic <- d[d$acute == 0, ]
vars <- c("age", "meanbp", "hrt", "resp", "temp", "crea", "sod", "wblc", "num.co", "male")
fit <- glm(reformulate(vars, "y"), family = binomial(), data = acute)
p <- predict(fit, chronic, type = "response")
X <- chronic[, vars]
sink(file.path(HERE, "worked_example.txt"))
v <- localize.update(chronic$y, p, X, B = 999, seed = 20261107)
print(v)
print(localize.size(chronic$y, p, X, R = 400, seed = 20261106))
sink()
ref <- read.csv(file.path(HERE, "size_apps.csv"))
print(ref[ref$model == "M0" & ref$setting == "transported, all chronic", c("part", "size", "chance", "lo", "hi")])
