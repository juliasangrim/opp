package ru.nsu.fit.g19202.jul.task1;

import java.util.Comparator;

public class WordComparator implements Comparator<Word> {
    @Override
    public int compare(Word a, Word b) {
            int result = b.get_freq().compareTo(a.get_freq());
            if (result == 0)
                result = a.get_word().compareTo(b.get_word());
            return result;
    }
}
