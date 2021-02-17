package ru.nsu.fit.g19202.jul.task1;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

public class DefaultWriter implements StatWriter {
    @Override
    public void writeStat(List<Word> list, WordComparator comp, int total, Writer writer) throws IOException {
        try (FileWriter out = new FileWriter("output.csv")) {
            list.sort(comp);
            for (Word element : list) {
                String percent = Float.toString((float) element.get_freq() / total * 100);
               writer.append(element.get_word())
                       .append("\t")
                       .append(element.get_freq().toString())
                       .append("\t")
                       .append(percent)
                       .append('\n');
            }
        }
    }
}
