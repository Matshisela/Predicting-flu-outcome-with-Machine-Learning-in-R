#Loading R packages
library(outbreaks)
library(tidyverse)
library(plyr)
library(mice)
library(caret)
library(purrr)

#Pre-processing
fluH7N9_china_2013$age[which(fluH7N9_china_2013$age == "?")] <- NA
fluH7N9_china_2013_gather <- fluH7N9_china_2013 %>%
  mutate(case_id = paste("case", case_id, sep = "_"),
         age = as.numeric(age)) %>%
  gather(Group, Date, date_of_onset:date_of_outcome) %>%
  mutate(Group = as.factor(mapvalues(Group, from = c("date_of_onset", "date_of_hospitalisation", "date_of_outcome"), 
                                     to = c("date of onset", "date of hospitalisation", "date of outcome"))),
         province = mapvalues(province, from = c("Anhui", "Beijing", "Fujian", "Guangdong", "Hebei", "Henan", "Hunan", "Jiangxi", "Shandong", "Taiwan"), to = rep("Other", 10)))

#Third variable added
levels(fluH7N9_china_2013_gather$gender) <- c(levels(fluH7N9_china_2013_gather$gender), "unknown")
fluH7N9_china_2013_gather$gender[is.na(fluH7N9_china_2013_gather$gender)] <- "unknown"
head(fluH7N9_china_2013_gather)


#Funtion
my_theme <- function(base_size = 12, base_family = "sans"){
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
      axis.title = element_text(size = 14),
      panel.grid.major = element_line(color = "grey"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "aliceblue"),
      strip.background = element_rect(fill = "lightgrey", color = "grey", size = 1),
      strip.text = element_text(face = "bold", size = 12, color = "black"),
      legend.position = "bottom",
      legend.justification = "top", 
      legend.box = "horizontal",
      legend.box.background = element_rect(colour = "grey50"),
      legend.background = element_blank(),
      panel.border = element_rect(color = "grey", fill = NA, size = 0.5)
    )
}


#Data Visualization
ggplot(data = fluH7N9_china_2013_gather, aes(x = Date, y = age, fill = outcome)) +
  stat_density2d(aes(alpha = ..level..), geom = "polygon") +
  geom_jitter(aes(color = outcome, shape = gender), size = 1.5) +
  geom_rug(aes(color = outcome)) +
  scale_y_continuous(limits = c(0, 90)) +
  labs(
    fill = "Outcome",
    color = "Outcome",
    alpha = "Level",
    shape = "Gender",
    x = "Date in 2013",
    y = "Age",
    title = "2013 Influenza A H7N9 cases in China",
    subtitle = "Dataset from 'outbreaks' package (Kucharski et al. 2014)",
    caption = ""
  ) +
  facet_grid(Group ~ province) +
  my_theme() +
  scale_shape_manual(values = c(15, 16, 17)) +
  scale_color_brewer(palette="Set1", na.value = "grey50") +
  scale_fill_brewer(palette="Set1")


ggplot(data = fluH7N9_china_2013_gather, aes(x = Date, y = age, color = outcome)) +
  geom_point(aes(color = outcome, shape = gender), size = 1.5, alpha = 0.6) +
  geom_path(aes(group = case_id)) +
  facet_wrap( ~ province, ncol = 2) +
  my_theme() +
  scale_shape_manual(values = c(15, 16, 17)) +
  scale_color_brewer(palette="Set1", na.value = "grey50") +
  scale_fill_brewer(palette="Set1") +
  labs(
    color = "Outcome",
    shape = "Gender",
    x = "Date in 2013",
    y = "Age",
    title = "2013 Influenza A H7N9 cases in China",
    subtitle = "Dataset from 'outbreaks' package (Kucharski et al. 2014)",
    caption = "\nTime from onset of flu to outcome."
  )


#Features
dataset <- fluH7N9_china_2013 %>%
  mutate(hospital = as.factor(ifelse(is.na(date_of_hospitalisation), 0, 1)),
         gender_f = as.factor(ifelse(gender == "f", 1, 0)),
         province_Jiangsu = as.factor(ifelse(province == "Jiangsu", 1, 0)),
         province_Shanghai = as.factor(ifelse(province == "Shanghai", 1, 0)),
         province_Zhejiang = as.factor(ifelse(province == "Zhejiang", 1, 0)),
         province_other = as.factor(ifelse(province == "Zhejiang" | province == "Jiangsu" | province == "Shanghai", 0, 1)),
         days_onset_to_outcome = as.numeric(as.character(gsub(" days", "",
                                                              as.Date(as.character(date_of_outcome), format = "%Y-%m-%d") - 
                                                                as.Date(as.character(date_of_onset), format = "%Y-%m-%d")))),
         days_onset_to_hospital = as.numeric(as.character(gsub(" days", "",
                                                               as.Date(as.character(date_of_hospitalisation), format = "%Y-%m-%d") - 
                                                                 as.Date(as.character(date_of_onset), format = "%Y-%m-%d")))),
         age = age,
         early_onset = as.factor(ifelse(date_of_onset < summary(fluH7N9_china_2013$date_of_onset)[[3]], 1, 0)),
         early_outcome = as.factor(ifelse(date_of_outcome < summary(fluH7N9_china_2013$date_of_outcome)[[3]], 1, 0))) %>%
  subset(select = -c(2:4, 6, 8))
rownames(dataset) <- dataset$case_id
dataset[, -2] <- as.numeric(as.matrix(dataset[, -2]))
head(dataset)

#Summary
summary(dataset$outcome)

#Imputing missing values
md.pattern(dataset)
dataset_impute <- mice(data = dataset[, -2],  print = FALSE)
datasets_complete <- right_join(dataset[, c(1, 2)], 
                                complete(dataset_impute, "long"),
                                by = "case_id") %>%
                                    select(-.id)

#compare the distributions of the five different imputed datasets:
datasets_complete %>%
  gather(x, y, age:early_outcome) %>%
  ggplot(aes(x = y, fill = .imp, color = .imp)) +
  facet_wrap(~ x, ncol = 3, scales = "free") +
  geom_density(alpha = 0.4) +
  scale_fill_brewer(palette="Set1", na.value = "grey50") +
  scale_color_brewer(palette="Set1", na.value = "grey50") +
  my_theme()


#Test, train and validation data sets
train_index <- which(is.na(datasets_complete$outcome))
train_data <- datasets_complete[-train_index, ]
test_data  <- datasets_complete[train_index, -2]


#In order to  model each imputed dataset separately, The nest() and map() functions are used
set.seed(263)
val_data <- train_data %>%
  group_by(.imp) %>%
  nest() %>%
  mutate(val_index = map(data, ~ createDataPartition(.$outcome, p = 0.7, list = FALSE)),
         val_train_data = map2(data, val_index, ~ .x[.y, ]),
         val_test_data = map2(data, val_index, ~ .x[-.y, ]))


#Random Forest
model_function <- function(df) {
  caret::train(outcome ~ .,
               data = df,
               method = "rf",
               preProcess = c("scale", "center"),
               trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, verboseIter = FALSE))
}

#Nested tibble
set.seed(42)
val_data_model <- val_data %>%
  mutate(model = map(val_train_data, ~ model_function(.x)),
         predict = map2(model, val_test_data, ~ data.frame(prediction = predict(.x, .y[, -2]))),
         predict_prob = map2(model, val_test_data, ~ data.frame(outcome = .y[, 2],
                                                                prediction = predict(.x, .y[, -2], type = "prob"))),
         confusion_matrix = map2(val_test_data, predict, ~ confusionMatrix(.x$outcome, .y$prediction)),
         confusion_matrix_tbl = map(confusion_matrix, ~ as.tibble(.x$table)))

#Comparing accuracy of models: Confusion matrices for the 5 datasets
val_data_model %>%
  unnest(confusion_matrix_tbl) %>%
  ggplot(aes(x = Prediction, y = Reference, fill = n)) +
  facet_wrap(~ .imp, ncol = 5, scales = "free") +
  geom_tile() +
  my_theme()

#prediction probabilities for correct and wrong predictions
val_data_model %>%
  unnest(predict_prob) %>%
  gather(x, y, prediction.Death:prediction.Recover) %>%
  ggplot(aes(x = x, y = y, fill = outcome)) +
  facet_wrap(~ .imp, ncol = 5, scales = "free") +
  geom_boxplot() +
  scale_fill_brewer(palette="Set1", na.value = "grey50") +
  my_theme()
