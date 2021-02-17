package ru.nsu.fit.g19202.jul.task1;

import java.io.IOException;
import java.io.Reader;
import java.util.Map;

public interface StatBuilder {
    public int readWords(Reader reader, Map<String, Word> words) throws IOException;
}