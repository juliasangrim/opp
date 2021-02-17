package ru.nsu.fit.g19202.jul.task1;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class Tests {

    @Test
    public void testComparator() {
        WordComparator comparator = new WordComparator();
        Word a = new Word("Alice", 1);
        Word b = new Word("Aright", 5);
        Word c = new Word("Bright", 5);
        assertEquals(1, comparator.compare(a, b));
        assertEquals(-1, comparator.compare(c, a));
        assertEquals(1, comparator.compare(c, b));
    }

    @Test
    public void testWriter() throws IOException {
        StringWriter writer = new StringWriter();
        StatWriter statWriter = new DefaultWriter();
        ArrayList<Word> list = new ArrayList<Word>();
        statWriter.writeStat(list, null, 0, writer);
        assertEquals("", writer.toString());
        Word word = new Word("Example", 1);
        list.add(word);
        statWriter.writeStat(list, null, 1, writer);
        assertEquals("Example\t1\t100.0\n", writer.toString());
    }

    @Test
    public void testStatBuilder() throws IOException {
        DefaultStatBuilder builder = new DefaultStatBuilder();
        StringReader reader = new StringReader("Hello world my friends, hello");
        Map<String,Word> words = new HashMap<>();
        int total = builder.readWords(reader, words);
        assertEquals(5, total);
        assertEquals(2, words.get("hello").get_freq());
        assertEquals(1, words.get("world").get_freq());
        assertEquals(1, words.get("my").get_freq());
        assertEquals(1, words.get("friends").get_freq());
    }
}
