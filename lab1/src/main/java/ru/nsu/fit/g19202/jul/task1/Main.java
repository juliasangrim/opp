package ru.nsu.fit.g19202.jul.task1;

import java.io.*;
import java.util.*;


public class Main {
    public  static void main(String[] args) {
        Map<String, Word> words = new HashMap<>();
        int total = 0;
        try (FileReader in = new FileReader(args[0])) {
            StatBuilder builder = new DefaultStatBuilder();
            total = builder.readWords(in, words);
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        var list = new ArrayList<>(words.values());
        WordComparator wordComp = new WordComparator();
        StatWriter statWriter = new DefaultWriter();
        try(Writer writer = getWriter(args)) {
            statWriter.writeStat(list, wordComp, total, writer);
        }
        catch (IOException e){
           e.printStackTrace();
        }
    }
    private static Writer getWriter(String[] args) throws IOException {
        if (args.length == 2) {
            return new FileWriter(args[1]);
        }
        return new OutputStreamWriter(System.out);
    }
}
