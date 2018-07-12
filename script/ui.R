shinyUI(
  fluidPage(
    titlePanel("voriPBPK"),
    sidebarPanel(
      sliderInput("Kpmu", "Kpmu", min=0.1, max=10, value=2.94),
      checkboxInput("optimKpmu", "Optimized Kpmu (0.78)", FALSE),
      sliderInput("BP", "BP", min=0.1, max=10, value=1),
      checkboxInput("optimBP", "Optimized BP (1.21)", FALSE),
      sliderInput("Peff", "Peff x 10^-4", min=0.01, max=1, value=0.145),
      checkboxInput("optimPeff", "Optimized Peff (0.04 x 10^-4)", FALSE),
      sliderInput("Adult_ClGu", "Adult Cl_Gu (mL/min/kg)", min=0, max=1, value=0),
      checkboxInput("optimAdult_ClGu", "Calculated Adult Cl_Gu (0.07 mL/min/kg)", FALSE),
      sliderInput("Pediatric_ClGu", "Pediatric Cl_Gu (mL/min/kg)", min=0, max=1, value=0),
      checkboxInput("optimPediatric_ClGu", "Calculated Pediatric Cl_Gu (0.19 mL/min/kg)", FALSE)
      ),
    mainPanel(
        plotOutput("myplot", height="1000", width="1000")
    )
  )
)
      
  #                      sliderInput("Age","Age (y)", min(nhanesData$AGE_YR), max = max(nhanesData$AGE_YR), value=30, step=1),
  #                      selectInput(
  #                        "Gender", "Gender", c("Male", "Female")
  #                      ),
  #                      sliderInput("Weight","Weight (kg)", min(nhanesData$BW), max = max(nhanesData$BW), value=73, step=1),
  #                      sliderInput("Height","Height (m)", min(nhanesData$HT/100), max = max(nhanesData$HT/100), value=1.76, step=0.01)
  #     ),
  #     
  #     conditionalPanel(condition="input.tabselected==2",
  #                      textInput("nSubj","Number of Subjects", value=10),
  #                      sliderInput("AgePop", "Age (y)", min = min(nhanesData$AGE_YR), max = max(nhanesData$AGE_YR), value = c(5,70), step=1),
  #                      textInput("femalePerc","Female Percentage", value=50),
  #                      sliderInput("WeightPop", "Weight (kg)", min = min(nhanesData$BW), max = max(nhanesData$BW), value = c(10,90), step=1),
  #                      sliderInput("HeightPop", "Height (m)", min = min(nhanesData$HT/100), max = max(nhanesData$HT/100), value = c(1,2), step=0.01)
  #                      
  #     ),
  #     conditionalPanel(condition="input.tabselected==3",
  #                      selectInput(
  #                        "prediction", "prediction", c("Poulin & Theil", "Berezhkovskiy", "Schmitt", "Rodgers & Rowland", "PK-Sim Standard")
  #                        #"prediction", "prediction", c("Poulin & Theil", "Berezhkovskiy", "Rodgers & Rowland")
  #                      ),
  #                      # conditionalPanel(condition="input.prediction == 'Poulin & Theil'",# || input.prediction == 'Berezhkovskiy' || input.prediction == 'Rodgers & Rowland' || input.prediction == 'Schmitt'",
  #                      #                  textInput("ka","ka", value=0.849),
  #                      #                  textInput("logP","logP", value=2.56),
  #                      #                  textInput("pKa","pKa", value=1.76),
  #                      #                  textInput("fup","fup", value=0.42),
  #                      #                  selectInput("compType", "Type", c("Acid", "Base"))),
  #                      # conditionalPanel(condition="input.prediction == 'Schmitt'",
  #                      #                  textInput("ka","ka", value=0.849),
  #                      #                  textInput("logP","logP", value=2.56),
  #                      #                  textInput("pKa","pKa", value=1.76),
  #                      #                  textInput("fup","fup", value=0.42),
  #                      #                  checkboxInput("Acidic", "Acidic", FALSE)),
  #                      # conditionalPanel(condition="input.prediction == 'Rodgers & Rowland'",
  #                      #                  textInput("ka","ka", value=0.849),
  #                      #                  textInput("logP","logP", value=2.56),
  #                      #                  textInput("pKa","pKa", value=1.76),
  #                      #                  textInput("fup","fup", value=0.42),
  #                      #                  selectInput("type", "Type", c("Acid", "Base"))),
  #                      conditionalPanel(condition="input.prediction == 'PK-Sim Standard'",
  #                                       textInput("ka","ka", value=0.849),
  #                                       textInput("logP","logP", value=2.56),
  #                                       textInput("fup","fup", value=0.42)),
  #                      conditionalPanel(condition="input.prediction != 'PK-Sim Standard'",
  #                                       textInput("ka","ka", value=0.849),
  #                                       textInput("logP","logP", value=2.56),
  #                                       textInput("pKa","pKa", value=1.76),
  #                                       textInput("fup","fup", value=0.42),
  #                                       selectInput("compType", "Type", c("Acid", "Base")))
  #     ),
  #     conditionalPanel(condition="input.tabselected==4",
  #                      textInput("Dose","Dose (mg/kg)", value=4),
  #                      selectInput(
  #                        "ROI", "ROI", c("IV", "PO")
  #                      ),
  #                      textInput("II","II (h)", value=12),
  #                      textInput("ADDL","ADDL", value=13),
  #                      textInput("Rate","Rate (dose/h)", value=4)
  #     ),
  #     conditionalPanel(condition="input.tabselected==5",
  #                      selectInput("simType", "Simulation:", choices=c("Individual", "Population")),
  #                      textInput("end", "End Time (h)", value=12),
  #                      textInput("delta", "delta", value=0.1),
  #                      checkboxInput("SS", "SS", FALSE)
  #                      
  #     ),
  #     conditionalPanel(condition="input.tabselected==5")
  #   ),
  #   mainPanel(
  #     tabsetPanel(
  #       tabPanel("Individual", value=1, conditionalPanel(condition="input.choice==1"), tableOutput("mytable1")),
  #       tabPanel("Population", value=2, conditionalPanel(condition="input.choice==2"), verbatimTextOutput("mysummary")),
  #       tabPanel("Compound", value=3, conditionalPanel(condition="input.choice==3"), tableOutput("mytable2")), 
  #       tabPanel("Administration Protocol", value=4, conditionalPanel(condition="input.choice==4"), plotOutput("myplot2")),
  #       tabPanel("Simulation", value=5, conditionalPanel(condition="input.choice==5"), plotOutput("myplot1")),
  #       id = "tabselected"
  #     )
  #   )
  # ))