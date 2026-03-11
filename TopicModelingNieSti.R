options(stringsAsFactors = FALSE)
library(quanteda)
library(topicmodels)
library(gutenbergr)
library(tidyverse)
library(quanteda.textstats)


Stirner <- readLines("Documents/OneDrive/Dokumente/Informatik/gutenberg(34580)", encoding = "UTF-8")
Nietzsche <- gutenberg_download(1998)
lemma_data <- read.csv("Documents/OneDrive/Dokumente/Informatik/RStudio/baseform en", encoding = "UTF-8")
stopwords_extended <- readLines("Documents/OneDrive/Dokumente/Informatik/RStudio/stopwords en", encoding = "UTF-8")

stirner_paragraphs <- Stirner[577:13369]
stirner_paragraphs <- stirner_paragraphs[nchar(stirner_paragraphs) > 50 & !is.na(stirner_paragraphs)]

nietzsche_paragraphs <- Nietzsche$text[677:13907]
nietzsche_paragraphs <- nietzsche_paragraphs[nchar(nietzsche_paragraphs) > 50 & !is.na(nietzsche_paragraphs)]

textdata <- data.frame(
  text = c(stirner_paragraphs, nietzsche_paragraphs),
  author = c(rep("Stirner", length(stirner_paragraphs)),
             rep("Nietzsche", length(nietzsche_paragraphs))),
  stringsAsFactors = FALSE
)


stopwords_extended <- c(stopwords_extended, "thou", "thy", "thee", "thyself", "zarathustra", "ye", "ah")


corpus_obj <- corpus(textdata, text_field = "text")


textdata_tokens <- corpus_obj %>% 
  tokens(remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE) %>% 
  tokens_tolower() %>% 
  tokens_replace(lemma_data$inflected_form, lemma_data$lemma, valuetype = "fixed") %>% 
  tokens_remove(pattern = stopwords_extended, padding = T)


textdata_collocations <- textstat_collocations(textdata_tokens, min_count = 25)
textdata_collocations <- textdata_collocations[1:250, ]
textdata_collocations <- textdata_collocations[!is.na(textdata_collocations$collocation), ]
textdata_tokens <- tokens_compound(textdata_tokens, textdata_collocations)

DTM <- textdata_tokens %>%
  tokens_remove("") %>%
  dfm() %>%
  dfm_trim(min_docfreq = 5)


top_terms_remove <- c("man", "thing", "one", "self")
DTM <- DTM[, !(colnames(DTM) %in% top_terms_remove)]


sel_idx <- rowSums(DTM) > 0
DTM <- DTM[sel_idx, ]
textdata <- textdata[sel_idx, ]


K <- 20
topicModel <- LDA(DTM,
                  K,
                  method = "Gibbs",
                  control = list(iter = 1000, seed = 1, verbose = 50, alpha = 0.1))

tmResult <- posterior(topicModel)
beta <- tmResult$terms    
theta <- tmResult$topics  


theta_df <- as.data.frame(theta)
theta_df$author <- docvars(DTM, "author")  
colnames(theta_df)[1:ncol(theta)] <- paste0("V", 1:ncol(theta))


top5termsPerTopic <- apply(beta, 1, function(x) {
  names(sort(x, decreasing = TRUE))[1:5]
})
topicNames <- apply(top5termsPerTopic, 2, paste, collapse = " ")

topic_mapping <- data.frame(
  topic = paste0("V", 1:K),
  topicName = topicNames,
  stringsAsFactors = FALSE
)


author_topic_means <- theta_df %>%
  group_by(author) %>%
  summarise(across(starts_with("V"), mean), .groups = "drop")

author_topic_long <- pivot_longer(author_topic_means,
                                  cols = -author,
                                  names_to = "topic",
                                  values_to = "mean") %>%
  left_join(topic_mapping, by = "topic") %>%
  filter(!is.na(topicName))  





ggplot(author_topic_long,
       aes(x = topicName, y = mean, fill = author)) +
  geom_col(stat = "identity", position = "dodge") +
  labs(
    x = "Topic (Top 5 Terms)",
    y = "Average Topic Proportion",
    title = "Topic Distribution: Nietzsche vs Stirner"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.title = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 40)
  )





library(LDAvis)
library("tsne")
svd_tsne <- function(x) tsne(svd(x)$u)
json <- createJSON(phi = beta, theta = theta, doc.length = rowSums(DTM),
                   vocab = colnames(DTM), term.frequency = colSums(DTM), mds.method = svd_tsne,
                   plot.opts = list(xlab = "", ylab = ""))
serVis(json)