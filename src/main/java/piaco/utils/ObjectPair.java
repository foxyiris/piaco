package piaco.utils;

import java.io.Serializable;

// Below page is a reference of this class.
// http://stackoverflow.com/questions/521171/a-java-collection-of-value-pairs-tuples
public class ObjectPair<T> implements Serializable
{

    private static final long serialVersionUID = 4438966336310741874L;

    private T f;
    private T s;

    public ObjectPair(T left, T right) {
        this.f = left;
        this.s = right;
    }

    public T getFirst(){
        return f;
    }
    
    public T getSecond() {
        return s;
    }

    @Override
    public int hashCode() {
        return f.hashCode() + s.hashCode();
    }

    @Override
    public boolean equals(Object o) {
        if (!(o instanceof ObjectPair)) return false;
        @SuppressWarnings("unchecked")
        ObjectPair<T> pairo = (ObjectPair<T>) o;
        boolean result = (this.f.equals(pairo.getFirst()) && this.s.equals(pairo.getSecond())) ||
                (this.f.equals(pairo.getSecond()) && this.s.equals(pairo.getFirst()));
                
        return result;
    }


    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        return sb.append("(").append(f.toString()).append(", ").append(s.toString()).append(")").toString();
    }
}
