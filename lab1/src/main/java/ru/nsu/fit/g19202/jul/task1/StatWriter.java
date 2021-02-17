package ru.nsu.fit.g19202.jul.task1;


import java.io.IOException;
import java.io.Writer;
import java.util.List;

public interface StatWriter {
    public void writeStat(List<Word> list, WordComparator comp, int total, Writer writer) throws IOException;
}