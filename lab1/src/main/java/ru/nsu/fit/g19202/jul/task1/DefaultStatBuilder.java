package ru.nsu.fit.g19202.jul.task1;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Map;

public class DefaultStatBuilder implements StatBuilder {
    public int readWords(Reader in, Map<String, Word> words) throws IOException {
        int total = 0;
            StringBuilder word = new StringBuilder();
            int c;
            while ((c = in.read()) != -1) {
                char cc = (char) c;
                if (Character.isLetterOrDigit(cc)) {
                    word.append(Character.toLowerCase(cc));
                }
                else if (word.length() > 0) {
                    String finalWord = word.toString();
                    words.compute(finalWord, (key, prev) -> {
                        if (prev == null) prev = new Word(finalWord, 1);
                        else  prev.change_freq();
                        return prev;
                    });
                    word.setLength(0);
                    total++;
                }
            }

            if (word.length() > 0) {
                String finalWord = word.toString();
                words.compute(finalWord, (key, prev) -> {
                    if (prev == null) prev = new Word(finalWord, 1);
                    else  prev.change_freq();
                    return prev;
                });
                word.setLength(0);
                total++;
            }
        return total;
    }
}
