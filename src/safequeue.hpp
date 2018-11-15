#include <mutex>
#include <queue>

template <class T>
class SafeQueue
{
private:
    std::queue<T> q;
    std::mutex    m;

public:

    inline void push(T t)
    {
        std::lock_guard<std::mutex> lock(m);
        q.push(t);
    }

    inline T pop()
    {
        std::lock_guard<std::mutex> lock(m);
        if (q.empty())
            return T();
        T val = q.front();
        q.pop();
        return val;
    }

    inline size_t size()
    {
        std::lock_guard< std::mutex > lock(m);
        return q.size();
    }

    inline bool empty()
    {
        std::lock_guard< std::mutex > lock(m);
        return q.empty();
    }
};
